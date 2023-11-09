### Setting
import re
import os	
import yaml
import pysam 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from rpy2.robjects import pandas2ri
from xml.etree.ElementTree import parse as parse_xml
from genomek.tools import chrom_cat_type_37 as chrom_cat_type

pandas2ri.activate()
%load_ext rpy2.ipython
%matplotlib inline 

pd.options.display.max_rows = 300
pd.options.display.max_columns = 100
pd.set_option('display.float_format', lambda x: '%.6f' % x)

! ls

df_sample = pd.read_csv("Sample_gem_name.csv", sep=',')
sample_ids = [str(x) for x in sorted(df_sample['Number'].unique())]
sample_config = {}
normal_fq_path = lambda sample_id, n: f"/home/users/euisoon123/DATA_GC/TBD220549_15816_20220824_WES_PBMC/Sample_sample_B{sample_id}/sample_B{sample_id}_{n}.fq.gz"
tumor_fq_path = lambda sample_id, n: f"/home/users/euisoon123/DATA_GC/TBD220549_16010_20220921_WES_Cancer/Sample_sample_T{sample_id}/sample_T{sample_id}_{n}.fq.gz"
bam_path = lambda sample_id: f"01_bam/{sample_id}/{sample_id}.splitmark.realigned.recal.bam"

for sample_id in sample_ids:
    sample_config[sample_id] = {}
    for type_id, fq_path in zip('NT', [normal_fq_path, tumor_fq_path]):
        sample_config[sample_id][type_id] = {}
        sample_config[sample_id][type_id]['id'] = type_id + sample_id
        sample_config[sample_id][type_id]['fq1'] = fq_path(sample_id, 1)
        sample_config[sample_id][type_id]['fq2'] = fq_path(sample_id, 2)
        sample_config[sample_id][type_id]['bam'] = bam_path(type_id + sample_id)

with open("sample_config.yaml", mode='w') as file:
    data = {}
    data['sample'] = sample_config
    yaml.dump(data, file, sort_keys=False, indent=4, default_style='"',)

with open('pon_bams.txt', 'w') as file:
    for bam_path in [type_dict['bam'] for sample_id, sample_dict in sample_config.items() for type_id, type_dict in sample_dict.items() if type_id == 'N']:
        print(os.path.abspath(bam_path), file=file)
        
with open('pon_keys.txt', 'w') as file:
    for type_id in [type_dict['id'] for sample_id, sample_dict in sample_config.items() for type_id, type_dict in sample_dict.items() if type_id == 'N']:
        print(type_id, file=file)

### Data Process
with open("sample_config.yaml", mode='r') as file:
    sample_config = yaml.safe_load(file)['sample']

import zipfile

def sampleId_to_fastqthroughput(sample_id):
    throughput_path = f"00_qc/fastq_throughput/{sample_id}/{sample_id}.tsv"
    fq1_throughput, fq2_throughput = pd.read_csv(f"00_qc/fastq_throughput/{sample_id}/{sample_id}.tsv", sep='\t', index_col=0, header=0).iloc[0,:]
    return fq1_throughput, fq2_throughput

def sampleId_to_fastqcStat(sample_id, fq, sample_yaml):
    fastqc_id = os.path.basename(sample_yaml[re.search(r"\d+", sample_id).group()][re.search(r"[NT]{1}", sample_id).group()][fq]).replace('.fq.gz', '_fastqc')
    with zipfile.ZipFile(f"00_qc/fastqc/{sample_id}/{fastqc_id}.zip") as thezip:
        with thezip.open(f'{fastqc_id}/fastqc_data.txt', mode='r') as thefile:
            for l in thefile:
                if l.startswith(b"Sequence length"):
                    length = l.decode('utf-8').split()[2]
                if l.startswith(b"%GC"):
                    GC = l.decode('utf-8').split()[1]
                    break
    return int(length), float(GC)

def sampleId_to_fastqscreenStat(sample_id, fq, sample_yaml):
    fastqc_id = os.path.basename(sample_yaml[re.search(r"\d+", sample_id).group()][re.search(r"[NT]{1}", sample_id).group()][fq]).replace('.fq.gz', '_screen')
    mouse_percent = pd.read_csv(f"00_qc/fastq_screen/{sample_id}/{fastqc_id}.txt", sep='\t', skiprows=1, header=0, index_col=0, usecols=['Genome', '%One_hit_one_genome', '%Multiple_hits_one_genome']).loc['Mouse',:].sum()
    return float(mouse_percent)

def sampleId_to_mosdepthInfo(sample_id):
    if os.path.exists(mosdepth_path := f"02_vcf/mosdepth_500/{sample_id}/{sample_id}.mosdepth.summary.txt") and os.path.getsize(mosdepth_path) > 0:
        mosdepth_summary = pd.read_csv(mosdepth_path, header=0, index_col=0, sep='\t')
        total_length, total_bases = mosdepth_summary.filter(like='region', axis=0).iloc[:-2,0:2].sum(axis=0)
        avg = total_bases/total_length
        try:
            mt_avg = float(mosdepth_summary.loc['MT', 'mean'])
        except:
            mt_avg = 0
        return avg, mt_avg
    else:
        return None, None

def sampleId_to_verifybamContam(sample_id):
    if os.path.exists(verifybam_path := f"00_qc/verifybamid2/{sample_id}.verifybam") and os.path.getsize(verifybam_path) > 0:
        try:
            return float(np.genfromtxt(verifybam_path, delimiter="\n", dtype=str)[-1].split(":")[1]) * 100
        except:
            return None
    else: 
        return None   
    
def sampleId_to_qprofilerResult(sample_id):
    if os.path.exists(qprofiler_path := f"00_qc/qprofiler/{sample_id}.qprofiler.xml") and os.path.getsize(qprofiler_path) > 0:
        qprofiler_xml = parse_xml(qprofiler_path).getroot().find("BAMReport")
        total_read = sum([int(x.attrib['count']) for x in qprofiler_xml.find("FLAG").find("ValueTally").findall("TallyItem")])
        dup_read = sum([int(x.attrib['count']) for x in qprofiler_xml.find("FLAG").find("ValueTally").findall("TallyItem") if 'd' in x.attrib['value']])
        unmap_read = sum([int(x.attrib['count']) for x in qprofiler_xml.find("FLAG").find("ValueTally").findall("TallyItem") if re.search(r"u|U{1}", x.attrib['value'])])
        qprofiler_isize = int(max([(x.attrib['count'], x.attrib['start']) for x in  qprofiler_xml.find("ISIZE").find("RG").find("RangeTally").findall("RangeTallyItem") if x.attrib['start'] != '0'], key=lambda x: int(x[0]))[1])
        qprofiler_isize_99 = xml_to_isize(qprofiler_xml, 0.99)
        qprofiler_dup_percent = dup_read/total_read * 100
        qprofiler_unmap_percent = unmap_read/total_read * 100
        return qprofiler_isize, qprofiler_isize_99, qprofiler_dup_percent, qprofiler_unmap_percent
    else:
        return None, None, None, None

def sampleId_to_samtoolsSingleton(sample_id, sample_yaml):
    if os.path.exists(singleton_path := sample_yaml[re.search(r"\d+", sample_id).group()][re.search(r"[NT]{1}", sample_id).group()]['bam'] + '.singleton') and os.path.getsize(singleton_path) > 0:
        return int(np.genfromtxt(singleton_path, dtype=int))
    else:
        return None
    
def sampleId_to_sequenzaCP(sample_id):
    def get_sequenzaCP(sequenza_cp_path):
        try:
            cellularity, ploidy = pd.read_csv(sequenza_cp_path, sep='\t').iloc[1:2,0:2].squeeze()
            return float(cellularity), float(ploidy)
        except:
            return 0, 0
    sequenza_old_path = f"02_vcf/sequenza/{sample_id}/{sample_id}_confints_CP.txt"
    sequenza_pjh_path = f"02_vcf/sequenza_pjh/{sample_id}/{sample_id}/{sample_id}_confints_CP.txt"
    sequenza_new_path = f"02_vcf/sequenza_filter/{sample_id}/{sample_id}/{sample_id}_confints_CP.txt"
    return sorted([get_sequenzaCP(sequenza_old_path), get_sequenzaCP(sequenza_pjh_path), get_sequenzaCP(sequenza_new_path)], key=lambda x: x[0], reverse=True)[0]

def sampleId_to_sequenzaHetCP(sample_id):
    def get_sequenzaCP(sequenza_cp_path):
        try:
            cellularity, ploidy = pd.read_csv(sequenza_cp_path, sep='\t').iloc[1:2,0:2].squeeze()
            return float(cellularity), float(ploidy)
        except:
            return 0, 0
    sequenza_het_path = f"02_vcf/sequenza_het/{sample_id}/{sample_id}_confints_CP.txt"
    return sorted([get_sequenzaCP(sequenza_het_path)], key=lambda x: x[0], reverse=True)[0]

def xml_to_isize(xml, percentile):
    '''
    receives xml and returns isize value of 99 percentile
    '''
    y, x = zip(*[(int(x.attrib['count']), int(x.attrib['start'])) for x in  xml.find("ISIZE").find("RG").find("RangeTally").findall("RangeTallyItem") if int(x.attrib['start']) > 0])
    for i in range(len(y)):
        if sum(y[:i]) > sum(y) * percentile:
            break
    return x[i-1:i][0]

sample_array = np.array(list(sample_config.keys()))

%%R -i sample_array -o result_array
library(glue)
result_array <- rep(0, length(sample_array))
names(result_array) <- sample_array
for(sample_id in sample_array) {
    load(glue("02_vcf/sequenza_filter/{sample_id}/{sample_id}_sequenza_extract.RData"))
    adr <- get(glue("{sample_id}_sequenza_extract"))
    adr <- adr$avg.depth.ratio
    result_array[sample_id] <- adr
}
adr_dict = {k:v for k,v in zip(sample_array, result_array)}

%%time
for patient_id, sample_dict in sample_config.items():
    for type_id, type_dict in sample_dict.items():
        if type_id not in "NT": continue
        sample_id = type_dict['id']

        type_dict['patient'] = patient_id
        type_dict['mosdepth_path'] = f"02_vcf/mosdepth_500/{sample_id}/{sample_id}.regions.bed.gz"
        type_dict['mosdepth_summary_path'] = f"02_vcf/mosdepth_500/{sample_id}/{sample_id}.mosdepth.summary.txt"
        if type_id in 'T':
            type_dict['sequenza_segment_path'] = f'02_vcf/sequenza_filter/{patient_id}/{patient_id}/{patient_id}_segments.txt'
            type_dict['sequenza_rdata_path'] = f"02_vcf/sequenza_filter/{patient_id}/{patient_id}/{patient_id}_sequenza_extract.RData"
            type_dict['varscan_loh_path'] = f"02_vcf/varscan/{patient_id}/{patient_id}.varscan.LOH.vcf.gz"
        type_dict['fq1_throughput'], type_dict['fq2_throughput'] = sampleId_to_fastqthroughput(sample_id)
        type_dict['fq1_length'], type_dict['fq1_GC'] = sampleId_to_fastqcStat(sample_id, 'fq1', sample_config)
        type_dict['fq2_length'], type_dict['fq2_GC'] = sampleId_to_fastqcStat(sample_id, 'fq2', sample_config)
        type_dict['fq1_mouse_onehit_percent'] = sampleId_to_fastqscreenStat(sample_id, 'fq1', sample_config)
        type_dict['fq2_mouse_onehit_percent'] = sampleId_to_fastqscreenStat(sample_id, 'fq2', sample_config)
        type_dict['avg_depth_mosdepth'], type_dict['MT_depth_mosdepth'] = sampleId_to_mosdepthInfo(sample_id)
        type_dict['cross_individual_contam'] = sampleId_to_verifybamContam(sample_id)
        type_dict['qprofiler_isize'], type_dict['qprofiler_isize_99'], type_dict['qprofiler_dup_percent'], type_dict['qprofiler_unmap_percent'] = sampleId_to_qprofilerResult(sample_id)
        type_dict['samflag_singleton'] = sampleId_to_samtoolsSingleton(sample_id, sample_config)
        if type_id in 'T':
            type_dict['sequenza_cellularity'], type_dict['sequenza_ploidy'] = sampleId_to_sequenzaCP(patient_id)
            type_dict['sequenza_avg_depth_ratio'] = float(adr_dict[patient_id])

%%time
from statsmodels.nonparametric.smoothers_lowess import lowess

for patient_id, sample_dict in sample_config.items():
    for key, value in sample_dict.items():
        if key not in "NT": continue
        df = pd.read_csv(value['mosdepth_path'], sep='\t', usecols=[0,3,4], names=['Chromosome', 'GC', 'cov'], dtype={'Chromosome': chrom_cat_type, 'GC': float, 'cov': float}).query("Chromosome != 'MT' and 0.3 <= GC <= 0.7")
        total_length, total_bases = pd.read_csv(value['mosdepth_summary_path'], header=0, index_col=0, sep='\t').filter(like='region', axis=0).iloc[:-2,0:2].sum(axis=0)
        avg = total_bases/total_length
        df['cov_normalized'] = df['cov']/avg
        x = np.arange(0.3, 0.7, 0.01)
        y = lowess(endog=df['cov_normalized'], exog=df['GC'], xvals=x)
        slope = np.diff(y)/np.diff(x)
        value['qc_gc_bias_variance'] = float(np.var(slope))
        value['qc_gc_bias_range'] = float(np.max(slope) - np.min(slope))

with open("total_meta_raw.yaml", 'w') as config:
    data = {}
    data['sample'] = sample_config
    yaml.dump(data, config, sort_keys=False, indent=4, default_style='"')

with open("total_meta_raw.yaml", 'r') as config:
    sample_yaml = yaml.safe_load(config)['sample']

meta_dict = {data_dict['id']: {k:v for k,v in data_dict.items() if 'path' not in k and k not in ['fq1', 'fq2', 'bam']} for group, types_dict in sample_yaml.items() for types, data_dict in types_dict.items() if types in 'NT'}
meta = pd.DataFrame.from_dict(meta_dict, orient='index')
meta

contam_df = pd.DataFrame.from_dict({patient_id:{type_id:type_dict['cross_individual_contam'] for type_id, type_dict in sample_dict.items() if type_id in "NT"} for patient_id, sample_dict in sample_yaml.items()})
contam_df = contam_df.T.sort_values(by=['N', 'T'], key=lambda x: pd.isna(x))[['N', 'T']].T
fig, ax = plt.subplots(figsize=(24,4))
sns.heatmap(contam_df, vmin=0, vmax=100, square=True, annot=True, fmt=".2f", ax=ax, cmap=sns.color_palette("rocket_r", as_cmap=True))

 

cellularity_df = pd.DataFrame.from_dict({patient_id:{type_id:type_dict['sequenza_cellularity'] for type_id, type_dict in sample_dict.items() if type_id in "T"} for patient_id, sample_dict in sample_yaml.items()})
cellularity_df = cellularity_df.loc[['T'], contam_df.columns]
fig, ax = plt.subplots(figsize=(24,4))
sns.heatmap(cellularity_df, vmin=0, vmax=1, annot=True, square=True, ax=ax, cmap=sns.light_palette("seagreen", as_cmap=True))

 

### NCM

import re
import genomek
from rpy2.robjects import pandas2ri
pandas2ri.activate()
%load_ext rpy2.ipython

%%R
draw_ncm <- function(matrix_path){
    output_corr_matrix <- read.delim(matrix_path)
    data = output_corr_matrix
    d3 <- as.dist((1 - data[,-1]))
    clust3 <- hclust(d3, method = "average")
    op = par(bg = "gray85")
    par(plt=c(0.1, 1, 0.2, 0.9))
    plot(clust3, lwd = 2, lty = 1,cex=0.8, xlab="Samples", sub = "",  ylab="Germline Distance (Pearson Correlation)",hang = -1, axes = FALSE)
    axis(side = 2, at = seq(0, 1, 0.2), labels = FALSE, lwd = 2)
    mtext(seq(0, 1, 0.2), side = 2, at = seq(0, 1, 0.2), line = 1,   las = 2)
}

with open("sample_config.yaml", mode='r') as file:
    sample_config = yaml.safe_load(file)['sample']

sample_list = {type_dict['id']:0 for sample, sample_dict in sample_config.items() for type_id, type_dict in sample_dict.items() if type_id in 'NT'}.keys()

%%time
! mkdir -p test/ngscheckmate/
sample_list = [x for x in sample_list]
genomek.qc.run_ncm([f"00_qc/ngscheckmate/vcf/{x}.vcf" for x in sample_list], outdir='test/ngscheckmate/')

%%R -w 4000 -h 500 -r 300 -p 6
draw_ncm("test/ngscheckmate/output_corr_matrix.txt")

 



