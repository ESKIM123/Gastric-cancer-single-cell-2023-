### Setting
import os
import re
import yaml
import pysam 
import numpy as np
import pandas as pd
import pyranges as pr
import matplotlib.pyplot as plt
import seaborn as sns
from rpy2.robjects import pandas2ri
from IPython.display import display, HTML
import genomek
from genomek.tools import chromosomes_37 as chromosomes
from genomek.tools import chrom_cat_type_37 as chrom_cat_type
from genomek.tools import chrom_sort_dict_37 as chrom_sort_dict
from genomek.visualization._basic_plots import genome_figure
from genomek.visualization._basic_plots import draw_scatter
from genomek.visualization._basic_plots import draw_bed
from genomek.visualization._basic_plots import draw_bnd
from genomek.visualization._basic_plots import draw_sv_on_fig


fasta = pysam.FastaFile("/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta")
sort_dict = {'N':0, 'T':1, 'X':2, 'O':3, 'NE':4, 'TE':5 ,'XE':6, 'Ne':7, 'Te':8 ,'Xe':9}
base_dict = {'A':'#35ab35', 'C':'#0000ff', 'G':'#d7862a', 'T':'#ff0000'}
hue_order  = ['C>T', 'C>A', 'C>G', 'T>C', 'T>A', 'T>G', 'others'] #,'G>A','G>C','G>T','T>A','T>C','T>G', 'not_snp']
hue_dict = {'C>A': [72/256,69/256,156/256], 'C>G': [1/256,1/256,1/256], 'C>T': [228/256,41/256,38/256], 
            'T>A': [168/256,53/256,144/256], 'T>C': [232/256,229/256,56/256], 'T>G': [110/256,173/256,43/256], 'others': 'grey'}
unimask_bed = pr.read_bed("/home/users/kimin/projects/00_Reference/unimask/um75-hs37d5.bed.gz")
hla_coord = ['6:29733420-29927046']
gene_bed = pr.read_bed("/home/users/pjh/References/igv_genomes/hg19/data/geneLocations_hg19.bed.gz").as_df().set_index("Name")
gene_bed['Chromosome'] = gene_bed['Chromosome'].str.replace('chr', '')

# Ignore future warnings from pandas
%matplotlib inline 
%load_ext autoreload
%autoreload 2

# R setting
pandas2ri.activate()
%load_ext rpy2.ipython

pd.options.display.max_rows = 300
pd.options.display.max_columns = 500
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 100

%%R
library(sequenza)
library(stringr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(circlize)
library(maftools)
library(glue)
library(GenomicRanges)
library(EnrichedHeatmap)

source("/home/users/kimin/projects/scripts/circos_kimin.R")

skip_variant = ['intron_variant', 'intergenic_variant', 'synonymous_variant']

oncogenes = pd.read_csv("~kimin/projects/00_Reference/cosmic_census/oncogenes_0506.txt", sep='\t', header=None)[0]
selected_genes = pd.read_csv("~kimin/projects/07_Gastric_Cancer/07_scripts/selected_genes.txt", sep='\t', header=None)[0]
rho_genes = np.array(['ARHGAP1','ARHGAP4','ARHGAP5','ARHGAP6','ARHGEF1','ARHGEF11','ARHGEF5','BAIAP2','DIAPH1','GSN','LIMK1','MYL2','OPHN1','PFN1','PIP5K1A','PIP5K1B','PPP1R12B','RHOA','ROCK1','TLN1','VCL'])
driver_genes = np.array(['TP53', 'RNF43', 'KRAS', 'PIK3CA', 'APC', 'MSH3', 'MSH6', 'MLH3', 'MLH1', 'POLD1', 'POLD4', 'ERBB2', 'ERBB3', 'ERBB4'])
TCGA_genes = np.array(['TP53','CDH1','SMAD4','PIK3CA','RHOA','ARID1A','KRAS','MUC6','APC','BCOR','EYA4','BNC2','RNF43','ABCA10','CTNNB1','MACF1','SMAD2','SOHLH2','RASA1','FAM46D','PLB1','CNGA4','EIF2C4','ERBB2','PTPRC'])
MSI_genes = np.array(['MSH3', 'PMS1', 'MLH3', 'EXO1', 'POLD1', 'POLD3', 'RFC1', 'RFC2', 'RFC3', 'RFC4', 'RFC5', 'PCNA', 'LIG1', 'RPA1', 'RPA2', 'RPA3', 'POLD2', 'POLD4', 'MLH1', 'MSH2', 'MSH6', 'PMS2'])

with open("total_meta_raw.yaml", 'r') as config:
    sample_yaml = yaml.safe_load(config)['sample']

meta_dict = {data_dict['id']: {k:v for k,v in data_dict.items() if 'path' not in k and k not in ['fq1', 'fq2', 'bam']} for group, types_dict in sample_yaml.items() for types, data_dict in types_dict.items() if types in 'NT'}
meta = pd.DataFrame.from_dict(meta_dict, orient='index')

def igv_load_group(group):
    ! echo new | nc localhost 60151
    for bam_path in dict(sorted({k:v['bam'] for k, v in sample_yaml[group].items() if k in 'NETEXEO'}.items(), key=lambda x: sort_dict[x[0]])).values():
        ! echo load {os.path.abspath(bam_path)} | nc localhost 60151
    ! echo collapse | nc localhost 60151

def igv_goto(chrom, pos):
    window_size = 150
    ! echo goto {str(chrom) + ":" + str(int(pos) - int(window_size/2)) + "-" + str(int(pos) + int(window_size/2))} | nc localhost 60151
    ! echo sort base | nc localhost 60151
#     ! echo collapse | nc localhost 60151
#     ! echo colorBy FIRST_OF_PAIR_STRAND | nc localhost 60151


def igv_load_bam(bam_path):
    ! echo load {os.path.abspath(bam_path)} | nc localhost 60151
    ! echo viewaspairs false | nc localhost 60151
    ! echo collapse | nc localhost 60151

def igv_goto_sv(chr1, pos1, chr2, pos2, svtype, window_size=2000):
    ! echo goto {str(chr1) + ":" + str(int(pos1) - int(window_size/2)) + "-" + str(int(pos1) + int(window_size/2))} {str(chr2) + ":" + str(int(pos2) - int(window_size/2)) + "-" + str(int(pos2) + int(window_size/2))}| nc localhost 60151
    if svtype == 'TRA':
        ! echo group mate_chromosome | nc localhost 60151
        ! echo sort matechr | nc localhost 60151
    else:
        ! echo group pair_orientation | nc localhost 60151
        ! echo viewaspairs | nc localhost 60151
        ! echo sort insertsize | nc localhost 60151

def igv_goto_seg(chrom, pos1, pos2):
    ! echo goto {str(chrom) + ":" + str(pos1) + "-" + str(pos2)} | nc localhost 60151
    ! echo sort base | nc localhost 60151

# ! echo new | nc localhost 60151
! echo squish | nc localhost 60151
# ! echo viewaspairs true | nc localhost 60151

### Filter shorts
%%time
dfs_short = {}
for pair in sample_yaml:
    print(pair)
    short_path = f"02_vcf/group_total/{pair}/{pair}.total.mle.ftr"
    df_short = pd.read_feather(short_path)
    df_short.loc[:, df_short.filter(regex="^(CALL|MUTECT|STRELKA|VARSCAN)_", axis=1).columns] = df_short.filter(regex="^(CALL|MUTECT|STRELKA|VARSCAN)_", axis=1).fillna(False)
    df_short.loc[:, df_short.filter(regex=".*_read_.*", axis=1).columns] = df_short.filter(regex=".*_read_.*", axis=1).fillna(0)
    for k in ['T']:
        df_short[f"vaf_highqual:::{k}"] = df_short[f"var_read_highqual:::{k}"] / df_short.filter(regex=f"^(ref)*(var)*(other)*_read_highqual:::{k}", axis=1).sum(axis=1).to_numpy()
    mask_bool = genomek.snv.variant_in_bed(df_short, pr.read_bed("/home/users/kimin/projects/00_Reference/unimask_kimin/hs37d5.kmask.bed")) 
    df_short['Mask_Okay'] = ~mask_bool
    mask_bool = genomek.snv.variant_in_bed(df_short, pr.read_bed("/home/users/kimin/projects/00_Reference/unimask/um75-hs37d5.bed.gz")) 
    df_short['Unimask_Okay'] = ~mask_bool
    df_short[f"Superpass"] = ( (df_short.filter(regex=f"^(MUTECT|STRELKA|VARSCAN)$").sum(axis=1)>=2) & \
                                   ( (df_short[f"var_read_all:::T"] - df_short[f"var_read_highqual:::T"] <= 10) | \
                                     (df_short[f"var_read_highqual:::T"]/df_short[f"var_read_all:::T"] >= 0.6) ) ) 
    df_short['mle_sum'] = df_short[['mle_het', 'mle_hom_ref', 'mle_hom_alt']].sum(axis=1)
    df_short['PoN_Okay'] = ((df_short['VEP_Existing_variation'] != '-') & (df_short['mle1'] > -200) & (df_short['var_read_all:::N'] <= 2 )) | \
                         ((df_short['VEP_Existing_variation'] == '-') & (df_short['pon_GC_vaf'] <= 0.1) & (df_short['var_read_all:::N'] <= 4 )) | \
                         ((df_short['VEP_Existing_variation'] == '-') & (df_short['pon_GC_vaf'] > 0.1) & (df_short['mle1'] > -200) & (df_short['var_read_all:::N'] <= 2 ))
    df_short['PoN_Okay2'] = ((df_short['VEP_Existing_variation'] != '-') & (df_short['pon_GC_vaf'] <= 0.01) & (df_short['var_read_all:::N'] <= 2 )) | \
                          ((df_short['VEP_Existing_variation'] != '-') & (df_short['pon_GC_vaf'] > 0.01) & (df_short['mle_sum'] > -10000) & (df_short['mle1'] > -300) & (df_short['var_read_all:::N'] <= 2 )) | \
                          ((df_short['VEP_Existing_variation'] == '-') & (df_short['pon_GC_vaf'] <= 0.01) & (df_short['var_read_all:::N'] <= 4 )) | \
                          ((df_short['VEP_Existing_variation'] == '-') & (df_short['pon_GC_vaf'] > 0.01) & (df_short['mle_sum'] > -10000) & (df_short['mle1'] > -300) & (df_short['var_read_all:::N'] <= 2 ))
    df_short['VarRead_Okay'] = (df_short.filter(regex="^var_read_highqual.*[TXOE]{1,2}$").sum(axis=1) / df_short.filter(regex="^var_read_all.*[TXOE]{1,2}$").sum(axis=1)) >= 0.7

    dfs_short[pair] = df_short
    print(df_short['VEP_VARIANT_CLASS'].value_counts())
    print(df_short['Superpass'].sum())

for sample, df_short in dfs_short.items():
    filter_bool = (df_short['pon_GC_vaf'] <= 0.01) | (df_short['var_read_highqual:::N'] <= 1)
    filter_bool = filter_bool & ((df_short['Unimask_Okay'] & (df_short['VEP_VARIANT_CLASS'] == 'SNV')) | (df_short['Mask_Okay'] & ((df_short['VEP_VARIANT_CLASS'] == 'deletion') | (df_short['VEP_VARIANT_CLASS'] == 'insertion'))))
    df_short['Filter'] = filter_bool
    df_short['Sample'] = sample
    print(sample, len(df_short), df_short['Filter'].sum())
    sns.distplot(df_short['vaf_highqual:::T'])
    plt.show()

df_short_total = pd.concat([df[df['Filter']] for df in dfs_short.values()], axis=0)

#### IGV Box
df_igv = df_short_total.sample(1)
sample = df_igv.iloc[0]['Sample']
igv_load_group(sample)
display(HTML(df_igv.to_html()))
chrom, pos = df_igv.iloc[0][['CHROM', 'POS']]
igv_goto(chrom, pos)
df_igv.filter(regex="Superpass")

### TMB
dfs_tmb = {}
for group, df in dfs_short.items():
    data = df.query('Filter')['VEP_VARIANT_CLASS'].value_counts()
    dfs_tmb[group] = {'SNV': data['SNV'], 
                      'INDEL': data['insertion']+data['deletion']}
df_tmb = pd.DataFrame.from_dict(dfs_tmb, orient='index').fillna(0)
sample_index = df_tmb.eval("Total = SNV + INDEL").sort_values(["Total", "INDEL"], ascending=False).index
print(sample_index)
df_tmb = df_tmb.reindex(sample_index, axis=0)
df_tmb.plot(kind='bar', stacked=True, figsize = (12,6), fontsize= 12)
plt.suptitle("Union Set")

 

### MAF
sample_index

gene_list = np.array(list(set(selected_genes).union(set(oncogenes)).union(set(TCGA_genes)).union(set(driver_genes)).union(set(MSI_genes))))

%%R
my_vs = c(
  "Frame_Shift_Del",
  "Frame_Shift_Ins",
  "Splice_Site",
  "Translation_Start_Site",
  "Nonsense_Mutation",
  "Nonstop_Mutation",
  "In_Frame_Del",
  "In_Frame_Ins",
  "Missense_Mutation"
)
flag_genes <- maftools:::flags(top = 20)
flag_genes <- c(flag_genes, 'POU6F2')
maf = read.maf(maf="03_analysis/maf/total.1231.maf", vc_nonSyn=my_vs)
maf <- filterMaf(maf=maf, genes=flag_genes)

sample_order = np.array(sample_index)

%%R -w 4 -h 4 -u in -r 300 -p 4 -i driver_genes -i sample_order
options(repr.plot.res=300)
oncoplot(maf=maf, top=50, genes=driver_genes, sampleOrder=sample_order,
         showTumorSampleBarcodes=TRUE, drawColBar=TRUE, removeNonMutated=FALSE, draw_titv=TRUE, keepGeneOrder=TRUE)

 


