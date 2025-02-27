#~/miniconda3/bin/python3

# ------- MATE/ Institute of Genetics and Biotechnology   ------- # 
# ------- Genomics and Bioinformatics Lab, M.Najjar, 2025 ------- #

## This script checks for the variants in the 10Kb upstream and Downstream a gene region and reports the complete set of variants at this surrounding region.
## Next, based on the case defined (H_L,L_H, H_M, L_M, or M_M), the script reports the matched variants genotype with the expression pattern at the given case 

import os,sys
import pandas as pd
import numpy as np
from argparse import ArgumentParser

parser = ArgumentParser(description="checks for the variants in the 10Kb upstream and Downstream, reports the matched variants genotype with the expression pattern at the given case" + "----- Genomics and Bioinformatics Lab, M.Najjar, 2025")

parser.add_argument('--input_file','-i',help="The input VCF file")
parser.add_argument('--case','-c',help="What is the expression CASE i.e. (H_L, H_M, L_H, L_M, or M_M)")
parser.add_argument('--annot_file','-a',help="The GTF annotation file")
parser.add_argument('--symbolized_file','-s',help="Case Symbolized file")

args = parser.parse_args()

CASE= args.case
LOG_FILE = args.symbolized_file
GTF_FILE = args.annot_file
VCF_FILE = args.input_file

Log2FC_DF = pd.read_csv(LOG_FILE, sep='\t')

## import gtf file

columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

# Read the GTF file 
GTF_109 = pd.read_csv(GTF_FILE, sep='\t', comment='#', header=None, names=columns)

# Extract the gene_id from the  (attribute column)
GTF_109.iloc[:, -1] =GTF_109.iloc[:, -1].str.extract(r'gene_id "([^"]+)"')


# Renaming the last column to 'gene_id'
GTF_109.rename(columns={GTF_109.columns[-1]: 'gene_id'}, inplace=True)
gtf_df = GTF_109.copy()



transcripts_df = gtf_df[gtf_df['feature'] == 'transcript']

# Group by gene_id and find the minimum  coordinates for each gene
gene_coords_df = transcripts_df.groupby('gene_id').agg({'seqname':'first','start': 'min', 'end': 'max'}).reset_index()


vcf_file = pd.read_csv(VCF_FILE, comment="#", header=None, sep="\t")
columns_names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO','INFO_GT','mother', 'child_1', 'child_2', 'child_3', 'child_4', 'child_5', 'child_6', 'child_7', 'child_8', 'father']


vcf_file.columns = columns_names
df = Log2FC_DF.copy()
FILTERED_DF = pd.DataFrame()
for index, row in df.iterrows():
    gene_name = row['Gene'].split('_')[0]
    filtered_row = gene_coords_df[(gene_coords_df['gene_id'] == gene_name)]
    FILTERED_DF = pd.concat([FILTERED_DF, filtered_row], ignore_index=True)

# Extracting the genotype information from the ""INFO_GT"
genotype_df = pd.DataFrame()

#Adding the chromosome and position columns
genotype_df["CHROM"] = vcf_file["CHROM"]
genotype_df["POS"] = vcf_file["POS"]

# Adding the individual names as columns
individuals = vcf_file.columns[9:]  # Assuming individual columns start from index 9
for ind in individuals:
    genotype_df[f"{ind}"] = vcf_file[ind].apply(lambda x: x.split(":")[0])


FILTERED_DF['start'] = pd.to_numeric(FILTERED_DF['start'])
FILTERED_DF['end'] = pd.to_numeric(FILTERED_DF['end'])
### dEfine the flanking

flanking_region_size = 10000


###  an empty list to store the filtered rows
filtered_rows = []

# Iteration  over rows in df
for _, row2 in FILTERED_DF.iterrows():

    gene_start = row2['start']
    gene_end = row2['end']
    gene_chromo = row2['seqname']

# Calculating flanking region's positions
    flank_start = gene_start - flanking_region_size
    flank_end = gene_end + flanking_region_size


    # Filter rows that fit within
    variants_within_range = genotype_df[
    (
        (genotype_df['POS'] >= flank_start) &
        (genotype_df['POS'] <= gene_start) &
        (genotype_df['CHROM'] == gene_chromo)
    ) |
    (
        (genotype_df['POS'] >= gene_end) &
        (genotype_df['POS'] <= flank_end) &
        (genotype_df['CHROM'] == gene_chromo)
    )
].copy()
    # # Adding additional columns to the filtered variants
    variants_within_range['gene_id'] = row2['gene_id']
    variants_within_range['start'] = row2['start']
    variants_within_range['end'] = row2['end']
    
    # ## To  filtered rows to the list
    filtered_rows.append(variants_within_range)

# Concatenating  the filtered rows ad  a single DataFrame
result_df = pd.concat(filtered_rows)
# Define a function to map genotype codes to alleles
def map_genotype(genotype):
    alleles = {'0': 'A', '1': 'B', '2': 'C'}  # Default alleles are for 0, 1, 2
    separator = '/' if '/' in genotype else '|'  ###  Determine a separator
    for i, allele in enumerate(genotype.split(separator)):
        if allele not in alleles:
            alleles[allele] = chr(65 + i)  # If we have more alleles
    alleles_list = [alleles[allele] for allele in genotype.split(separator)]
    return ''.join(alleles_list)

# Apply the function to the genotype columns
result_df['mother'] = result_df['mother'].apply(map_genotype)
result_df['child_1'] = result_df['child_1'].apply(map_genotype)
result_df['child_2'] = result_df['child_2'].apply(map_genotype)
result_df['child_3'] = result_df['child_3'].apply(map_genotype)
result_df['child_4'] = result_df['child_4'].apply(map_genotype)
result_df['child_5'] = result_df['child_5'].apply(map_genotype)
result_df['child_6'] = result_df['child_6'].apply(map_genotype)
result_df['child_7'] = result_df['child_7'].apply(map_genotype)
result_df['child_8'] = result_df['child_8'].apply(map_genotype)
result_df['father'] = result_df['father'].apply(map_genotype)

result_df.reset_index(drop=True, inplace=True)
result_df.drop_duplicates(inplace=True)
columns_to_compare = ['mother', 'child_1', 'child_2', 'child_3', 'child_4', 'child_5', 'child_6', 'child_7', 'child_8', 'father']
all_same_values = result_df[columns_to_compare].eq(result_df[columns_to_compare].iloc[:, 0], axis=0).all(axis=1)
rows_with_same_values = result_df[all_same_values]

merged_df = pd.merge(result_df, rows_with_same_values, how='left', indicator=True)

###  Filtering out the rows that are only present
Filtered_DF = merged_df[merged_df['_merge'] == 'left_only'].drop(columns='_merge')
Filtered_DF.to_csv(CASE+"_complete.tsv", sep='\t', index=False)

## Matching Genotype with phenotype 
if CASE in ['H_L','L_H']:
    df = Filtered_DF.copy()
    filtered_df = df[(df['mother'].isin(['AA', 'BB'])) &
                 (df['father'].isin(['AA', 'BB'])) &
                 (df['child_1'] == 'AB') &
                 (df['child_2'] == 'AB') &
                 (df['child_3'] == 'AB') &
                 (df['child_4'] == 'AB') &
                 (df['child_5'] == 'AB') &
                 (df['child_6'] == 'AB') &
                 (df['child_7'] == 'AB') &
                 (df['child_8'] == 'AB')]
    filtered_df.to_csv(CASE+"_matched_pattern.tsv", sep='\t', index=False)

elif CASE in ['H_M','L_M','M_M']:
    final_df = pd.DataFrame()
    df2 = Log2FC_DF.copy()
    df1 = Filtered_DF.copy()
    POSSIBLE_HET = ['AB','AC','BC','CD','AD','BD']
    # Loop over rows of the second DataFrame
    for index, row in df2.iterrows():
        ### Extracting the patterns from  this current row
        mother_like_individuals = set()
        L_like_individuals = set()
        H_like_individuals = set()
        MUM = row.iloc[1:].iloc[0]
        mother_like_individuals.add('mother')
        for col in row.index[2:]:
            if row[col] == MUM:
                mother_like_individuals.add(col)
            if row[col] == 'H':
                L_like_individuals.add(col)
            if row[col] == 'L':
                L_like_individuals.add(col)
        mother_like_individuals_sorted = sorted(list(mother_like_individuals))
        filtered_df1 = df1[(df1['gene_id'] == row['Gene'])]
        filtered_df2 = pd.DataFrame(columns=filtered_df1.columns)
        for index, row in filtered_df1.iterrows():
            mother_like_individuals_GT = set()
            L_like_GT = set()
            H_like_GT = set()
            MUM = row.iloc[2:].iloc[0]
            DAD = row.iloc[11:].iloc[0]
            POS = row.iloc[1:].iloc[0]
            mother_like_individuals_GT.add('mother')
            for col in row.index[2:]:
                if row[col] == MUM:
                    mother_like_individuals_GT.add(col)
                if col in H_like_individuals:
                    H_like_GT.add(row[col])
                if col in L_like_individuals:
                    L_like_GT.add(row[col])
            mother_like_individuals_GT_sorted = sorted(list(mother_like_individuals_GT))
            if mother_like_individuals_GT_sorted == mother_like_individuals_sorted: #matching the genotype with phenotype baseon the likelihood with mother
                if CASE == 'M_M':
                    if L_like_GT and len(L_like_GT) !=1:
                        #print("skipped more than gt in L", POS) #only one genotype should be responsible for L or H in M_M case 
                        continue
                    elif H_like_GT and len(H_like_GT) !=1:
                        #print("skipped more than gt in H", POS)
                        continue
                    elif (DAD == MUM) and (DAD in POSSIBLE_HET):     ###insure the parents are the same genotype

                        row_df = pd.DataFrame(row).T
                    else:
                        print("cases is M_M and skipped",POS)
                        continue
                else:
                    if (DAD in POSSIBLE_HET) and (MUM in POSSIBLE_HET): # Both parents can not be heterozygous
                        print("case is not (H_M/L_M) skipped",POS)
                        continue
                    else:
                        #print("Added",POS,MUM,DAD,sep='\t')
                        row_df = pd.DataFrame(row).T
                filtered_df2 = pd.concat([filtered_df2, row_df], ignore_index=True)
        final_df = pd.concat([final_df, filtered_df2])
    final_df.reset_index(drop=True, inplace=True)
    final_df.to_csv(CASE+"_matched_pattern.tsv", sep='\t', index=False)

