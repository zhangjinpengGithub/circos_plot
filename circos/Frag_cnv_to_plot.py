import pandas as pd
import sys,os
import subprocess
import argparse
import pandas as pd
import numpy as np
import argparse
from PIL import Image
def run_command(cmd):
#       print(cmd)
        return_code = subprocess.call(cmd, shell=True)
def frag_plot(frag_name,script):
    output=frag_name.split("/")[-1].split(".")[0]+".frag"
    workscript=script
    frature=pd.read_csv(workscript+"/frag8.1",header=0)
    cancersample=pd.read_csv(frag_name,header=0)
    frag_bin=pd.read_table(workscript+"/frag.bin",header=0)
    cancersample['SampleID']='sample'
    cancers=pd.melt(cancersample,id_vars=['SampleID'])
    cancers.columns=['flag','bin','value']
    cancers['bin']=cancers['bin'].map(lambda x:str(x).replace("ratio","bin"))
    cancerfrag=pd.merge(cancers,frag_bin,on=['bin'])
    cancerfrag['SampleID']='sample'
    a2=cancerfrag[['SampleID','bin','value','seqnames','start','end','flag']]
    frature['SampleID']=frature['Type']
    frature['flag']=frature['Type']
    frature=frature[a2.columns.to_list()]
    a3=pd.concat([frature,a2])
    a3['flag']=a3['flag'].map(lambda x:str(x).replace("Cancer","cancer_base").replace("Healthy","healthy_base"))
    a3.to_csv(output_dir+"/"+output,sep='\t',index=False)
###RR
    cmd="/dssg/home/sheny/MyProject/Mercury2/conda_env/8e6ace282ebb2428b18c350af5216cda_/bin/Rscript "+workscript+"/frag_plot.R "+output_dir+"/"+output +" "+workscript+"/cnv.base"
    print(cmd)
    run_command(cmd)
#    run_command("/dssg/home/zhangjp/miniconda3/bin/python "+workscript+"/pdf2png.py " +output+".pdf "+output)
def cnvflag(x):
    cnv_region=x['index']
    if(x[0]>0 and x[0]>ci[cnv_region][1]):
        return 'Amp'
    elif(x[0]<0 and x[0]<ci[cnv_region][0]):
        return 'Del'
    else:
        return 'Normal'

def cnv_plot(i,script):
    workscript=script
    CNV_Healthy=pd.read_csv(workscript+"/Cnv.Healthy.base",index_col=0 )
    df=CNV_Healthy.dropna(axis=1)
    lower = np.percentile(df, q=0.5, axis=0)
    upper = np.percentile(df, q=99.5, axis=0)
    global ci
    ci = pd.DataFrame([lower, upper], index=['lower', 'upper'])
    ci.columns=df.columns
    sample=pd.read_csv(i,index_col=0)
    sample=sample.dropna(axis=1)
    sample=sample.mean().reset_index()
    sample[['cnv','chr','start','end']]=sample['index'].str.split('.', expand=True)
    sample['chr']=sample['chr'].map(lambda x:"hs"+str(x))
    sample_filter=sample.copy()
    sample=sample[['chr','start','end',0]]
    sample.columns=['chr','start','end','log2ratio']
    os.chdir(output_dir)
    sample.to_csv("./sample.cnv.txt", sep = '\t', index = False,header=None)
    sample_filter.loc[:,"flag"] = sample_filter.apply(cnvflag,axis=1)
    sample_filter=sample_filter.query("flag!='Normal'")
    diffs = {}
    for x in sample_filter.itertuples():
        mid=ci[x[1]].mean()
        diff=abs(x[2]-mid)
        diffs[x[1]]=diff
    d_sorted = sorted(diffs.items(), key=lambda x: x[1], reverse=True)[:15]
    index_list=[]
    for idx in d_sorted:
        index_list.append(idx[0])
    sample_filter=sample_filter.query("index in @index_list")
    sample_filter=sample_filter[['chr','start','end',0]]
    sample_filter.to_csv("./sample.cnv_sample_filter.txt", sep = '\t', index = False,header=None)
    ###
    run_command(f"cp {workscript}/karyotype.human_chr.txt {workscript}/Cnv_cancer.cnv.txt {workscript}/Cnv_Healthy.cnv.txt ./")
    run_command("/dssg/home/wutt/anaconda3/envs/zjptest/bin/circos  -conf  "+workscript+"/circos_sample_2023.conf ") #https://circos.ca/software/download/circos/ 
    run_command("/dssg/home/zhangjp//miniconda3/envs/NEXTFlow/bin/cairosvg circos.svg -o circos.png -s 3 ") #pip install CairoSVG
    im = Image.open("./circos.png")   ###from PIL import Image
    cropped_im = im.crop((2000, 2000, 7000, 7000))
    cropped_im.save(i.split("/")[-1].split(".")[0]+".cnv.png")
   # run_command("rm karyotype.human_chr.txt Cnv_cancer.cnv.txt Cnv_Healthy.cnv.txt circos.png circos.svg ")
def main(args):
    global output_dir
    cnv=args.cnv
    frag=args.frag
    script=args.script
    output_dir=args.output_dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if(os.path.exists(args.script)==False ):
        sys.exit(' Error!Please check script path')
    print(script)
    run_command("cp "+frag+" "+output_dir)
    run_command("cp "+cnv+" "+output_dir)
    frag_plot(frag,script)
    cnv_plot(cnv,script)
if __name__ == "__main__":
        parser = argparse.ArgumentParser(
                prog="python Frag_cnv_to_plot.py  --cnv SampleID.cnv.csv --frag SampleID_fragment_ratio.csv --script /dssg/home/zhangjp/tmp/L008_demo/script --output_dir Output_DIR ",
                formatter_class=argparse.RawTextHelpFormatter,
                description='''
Usage: python Frag_cnv_to_plot.py  [Options]
''')
        parser.add_argument("--cnv")
        parser.add_argument("--frag")
        parser.add_argument("--output_dir") 
        parser.add_argument("--script",default='/dssg/home/zhangjp/tmp/L008_demo/script')
        args = parser.parse_args()
        main(args)
        print("All tasks done")

