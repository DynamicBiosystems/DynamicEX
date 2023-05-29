import json
import argparse
import pandas as pd
import numpy as np
from collections import defaultdict
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()
parser.add_argument('--STARSummary',required=True)
parser.add_argument('--STARLog',required=True)
parser.add_argument('--RnaSeqMetrics',required=True)
parser.add_argument('--UMIpreCell',required=True)
parser.add_argument('--sampleinfo',required=True)
parser.add_argument("--report_temp",required=True)
parser.add_argument("--output",required=True)
parser.add_argument("--fastp",required=True)
parser.add_argument("--MetricsSummary",required=True)

args = parser.parse_args()

def trans(i):
    if i <= 100000:
        return str(i).replace("000","")+"k"
    else:
        return str(i).replace("000000","")+"M"
        
def Barcode_rank_plot(UMIpreCell,cells,output):
    sns.set(style="white")
    plt.figure(figsize=(7.3,6.5),dpi=200)
    df = pd.read_csv(UMIpreCell,header=None)
    cells = int(cells.replace(",",""))
    df['tag'] = ["Cells" if i < cells else "Backgroud" for i in range(len(df))]
    index = 0
    for tag in ['Cells',"Backgroud"]:
        dft = df.query(f"tag=='{tag}'")
        if tag == "Cells":
            plt.plot(range(len(dft)),dft[0].values,color="#153057",lw=2.5,label="Cells")
            index = len(dft)
        else:
            plt.plot(range(index,len(dft)+index),dft[0].values,color="grey",lw=2.5,label="Backgroud")
    plt.yscale('log')
    plt.xscale('log')
    plt.legend(edgecolor="white")
    #####坐标轴
    plt.gca().spines["right"].set_alpha(.0)
    plt.gca().spines["top"].set_alpha(.0)
    plt.title("Barcode Rank Plot",fontweight ="bold",fontsize=20)
    plt.xlabel("Barcodes",fontweight ="bold",fontsize=20)
    plt.ylabel("UMI counts",fontweight ="bold",fontsize=20)
    xticks = [10**i for i in range(0,100) if 10**i <len(df)]
    xtickslabel = [i if i <=100 else trans(i) for i in xticks]
    yticks = [10**i for i in range(0,100) if 10**i <df[0].max()]
    ytickslabel = [i if i <=100 else trans(i) for i in yticks]
    plt.xticks(xticks,xtickslabel)
    plt.yticks(yticks,ytickslabel)
    plt.grid()
    plt.tight_layout()
    src = output.parent/'src'
    src.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(src/"barcode_rank_plot.png"))
    
def cell_and_sequencing_summary(STARSummary):
    summary = {"cell":{},"sequencing":{}}
    df = pd.read_csv(STARSummary,index_col=0,header=None)
    df = df.dropna()
    df.columns=['index']
    d = df.to_dict()['index']
    summary['sequencing']["Reads Mapped to Genome"] = f"{d['Reads Mapped to Genome: Unique+Multiple']:.2%}"
    summary['sequencing']['Number of Reads'] = f"{int(d['Number of Reads']):,}"
    summary['sequencing']['Valid Barcodes'] = f"{d['Reads With Valid Barcodes']:.2%}"
    summary['sequencing']['Sequencing Saturation'] = f"{d['Sequencing Saturation']:.2%}"
    if d.get('Q30 Bases in CB+UMI'):
        summary['sequencing']['CB+UMI'] = f"""
<tr>
  <th>Q30 Bases in CB+UMI</th>
  <th>{d['Q30 Bases in CB+UMI']:.2%}</th>
</tr>
"""
    else:
        summary['sequencing']['CB+UMI'] = ""
    summary['sequencing']['Q30 Bases in RNA Read'] = f"{d['Q30 Bases in RNA read']:.2%}"
    summary['cell']['Estimated Number of Cells'] = f"{int(d['Estimated Number of Cells']):,}"
    summary['cell']['Fraction Reads in Cells'] = f"{d['Fraction of Unique Reads in Cells']:.2%}"
    summary['cell']['Mean Reads per Cell'] = f"{int(int(d['Number of Reads'])/int(d['Estimated Number of Cells'])):,}"
    summary['cell']['Median UMI Counts per Cell'] = f"{int(d['Median UMI per Cell']):,}"

    if 'Median Gene per Cell' in d:
        summary['cell']['Median Genes per Cell'] = f"{int(d['Median Gene per Cell']):,}"
        summary['cell']['Total Genes Detected'] = f"{int(d['Total Gene Detected']):,}"
    else:
        summary['cell']['Median Genes per Cell'] = f"{int(d['Median GeneFull per Cell']):,}"
        summary['cell']['Total Genes Detected'] = f"{int(d['Total GeneFull Detected']):,}"

    return summary

def mapping_summary(STARLog, RnaSeqMetrics):
    d={}
    summary = {}
    with open(STARLog, 'r',encoding='utf-8') as fh:
        for line in fh:
            if 'Number of input reads' in line:
                summary['Number of input reads'] = int(
                    line.strip().split('\t')[-1])
            if 'Uniquely mapped reads number' in line:
                summary['Uniquely mapped reads number'] = int(
                    line.strip().split('\t')[-1])
            if 'Number of reads mapped to multiple loci' in line:
                summary['Number of reads mapped to multiple loci'] = int(
                    line.strip().split('\t')[-1])
    with open(RnaSeqMetrics, 'r',encoding='utf-8') as fh:
        while True:
            line = fh.readline().strip()
            if line.startswith('total alignments'):
                summary['total alignments'] = int(line.split()[-1].replace(
                    ',', ''))
            if line.startswith('reads aligned'):
                summary['reads aligned'] = int(line.split()[-1].replace(
                    ',', ''))
            if line.startswith('aligned to genes'):
                summary['aligned to genes'] = int(line.split()[-1].replace(
                    ',', ''))
            if line.startswith('no feature assigned'):
                summary['no feature assigned'] = int(line.split()[-1].replace(
                    ',', ''))
            if line.startswith('exonic'):
                summary['exonic'] = int(line.split()[-2].replace(',', ''))
            if line.startswith('intronic'):
                summary['intronic'] = int(line.split()[-2].replace(',', ''))
            if line.startswith('intergenic'):
                summary['intergenic'] = int(line.split()[-2].replace(',', ''))
                break
    #d['Reads Mapped to Genome'] = f"{summary['reads aligned']/summary['Number of input reads']:.2%}"
    intergenic = summary['intergenic']/summary['Number of input reads']
    intronic = summary['intronic']/summary['Number of input reads']
    exonic = summary['exonic']/summary['Number of input reads']
    d['Reads Mapped Confidently to Genome'] = f"{intergenic+intronic+exonic:.2%}"
    d['Reads Mapped Confidently to Intergenic Regions'] = f"{intergenic:.2%}"
    d['Reads Mapped Confidently to Intronic Regions'] = f"{intronic:.2%}"
    d['Reads Mapped Confidently to Exonic Regions'] = f"{exonic:.2%}"
    return d

def Sample_summary(config):
    summary = {}
    d = json.load(open(config))
    summary['Sample ID'] = d['sample']
    summary['Sample Description'] = d.get("Sample Description","")
    summary['Chemistry'] = "Single Cell 3' Decoder1.0"
    summary['Include introns'] = str(d['include_introns'])
    summary['Transcriptome'] = Path(d['transcriptome']).parent.name
    summary['Pipeline Version'] = "DynamicEX 1.0.0"
    return summary

def fastp_qc(fastp):
    summary={}
    d = json.load(open(fastp))
    summary["gc_content"] = f"{d['summary']['before_filtering']['gc_content']:.2%}"
    summary["too_short_reads"] = f"{d['filtering_result']['too_short_reads']:,}"
    summary["too_many_N_reads"] = f"{d['filtering_result']['too_many_N_reads']:,}"
    summary['q20_rate'] = f"{d['summary']['before_filtering']['q20_rate']:.2%}"
    return summary
    
STARSummary = args.STARSummary
RnaSeqMetrics = args.RnaSeqMetrics
sampleinfo = args.sampleinfo
STARLog = args.STARLog
Temp = args.report_temp
UMIpreCell = args.UMIpreCell
output = Path(args.output)
fastp = args.fastp

total_summary = {}
total_summary["qc"] = fastp_qc(fastp)
total_summary.update(cell_and_sequencing_summary(STARSummary))
cells = total_summary['cell']['Estimated Number of Cells']
Barcode_rank_plot(UMIpreCell,cells,output)
total_summary['mapping'] = mapping_summary(STARLog, RnaSeqMetrics)

metrics_summary = []

ms = {}
for key in total_summary:
    ms.update(total_summary[key])
trans = {"GC Content":"gc_content","Q20 Bases in RNA Read":"q20_rate"}
print(ms)
for key in ["Number of Reads","Valid Barcodes","Sequencing Saturation","GC Content","Q20 Bases in RNA Read","Q30 Bases in RNA Read","Estimated Number of Cells","Fraction Reads in Cells","Mean Reads per Cell","Median UMI Counts per Cell","Median Genes per Cell","Total Genes Detected","Reads Mapped to Genome","Reads Mapped Confidently to Genome","Reads Mapped Confidently to Intergenic Regions","Reads Mapped Confidently to Intronic Regions","Reads Mapped Confidently to Exonic Regions"]:
    key2 = trans.get(key,key)
    metrics_summary.append([key,str(ms[key2]).replace(',','')])
print(pd.DataFrame(metrics_summary))
pd.DataFrame(metrics_summary).to_csv(args.MetricsSummary,index=False,header=False)

total_summary['sample'] = Sample_summary(sampleinfo)

f = open(Temp,encoding='utf-8').read()
for key in total_summary:
    for info in total_summary[key]:
        f = f.replace(f"$${info}$$",str(total_summary[key][info]))
        
with open(output,"w",encoding='utf-8') as out:
    out.write(f)
