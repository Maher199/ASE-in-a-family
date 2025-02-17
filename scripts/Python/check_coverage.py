import pysam
import pandas as pd
import sys

## This  script checks the exons in parents' bam files and reports the genes that have at least one exon that has a continuous coverage more than 5 bp 

### These arguments shouls be provided: 1- list of genes to check with "Gene" in the header
#					2- GTF file
#					3- Bam file for the mother sample
#					4- Bam file for the father sample

genes_to_check = sys.argv[1]
gtf_file = sys.argv[2]
bam_file_mother = sys.argv[3]
bam_file_mother = sys.argv[4]

## Prepare GTF_DF
columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
# Read the GTF file into a pandas DataFrame
GTF_109 = pd.read_csv(gtf_file, sep='\t', comment='#', header=None, names=columns)
# Extract the gene_id from the last column (attribute column)
GTF_109.iloc[:, -1] =GTF_109.iloc[:, -1].str.extract(r'gene_id "([^"]+)"')
# Rename the last column to 'gene_id'
GTF_109.rename(columns={GTF_109.columns[-1]: 'gene_id'}, inplace=True)

Genes_To_Check = pd.read_csv(genes_to_check,sep='\t')

FILTERED_GTF_DF = pd.DataFrame()
for index, row in Genes_To_Check.iterrows():
    gene_name = row['Gene'].split('_')[0]
    filtered_row = GTF_109[(GTF_109['gene_id'] == gene_name) & (GTF_109['feature'] == 'exon')]
    FILTERED_GTF_DF = pd.concat([FILTERED_GTF_DF, filtered_row], ignore_index=True)

Final_DF = pd.DataFrame(columns = ['Gene_ID', 'Exon_Start', 'Exon_End', 'Position', 'N_Count'])
for idx, row in FILTERED_GTF_DF.iterrows():
    EXON = []
    samfile = pysam.AlignmentFile(bam_file_mother, "rb" )
    CHROM = row['seqname']
    START = row['start']
    END = row['end']
    for pileupcolumn in samfile.pileup(CHROM, START, END):
        if pileupcolumn.pos >= START and pileupcolumn.pos <= END:
            new_row = [row['gene_id'],START,END, pileupcolumn.pos, pileupcolumn.n]
            EXON.append(new_row)
    samfile.close()
    EXON_df = pd.DataFrame(EXON)
    if not EXON_df.empty:
        EXON_df.columns = ['Gene_ID', 'Exon_Start', 'Exon_End','Position', 'N_Count']
        all_greater_than_5 = (EXON_df['N_Count'] > 5).all()
        if all_greater_than_5:
            Final_DF = pd.concat([Final_DF, EXON_df], ignore_index=True)
        else:
            print(row['gene_id'],START,END, "  Not found in Mother, Searching in Father")
            EXON = []
            samfile = pysam.AlignmentFile(bam_file_mother, "rb" )
            for pileupcolumn in samfile.pileup(CHROM, START, END):
                if pileupcolumn.pos >= START and pileupcolumn.pos <= END:
                    new_row = [row['gene_id'],START,END, pileupcolumn.pos, pileupcolumn.n]
                    EXON.append(new_row)
            samfile.close()
            EXON_df = pd.DataFrame(EXON)
            if not EXON_df.empty:
                EXON_df.columns = ['Gene_ID', 'Exon_Start', 'Exon_End','Position', 'N_Count']
                all_greater_than_5 = (EXON_df['N_Count'] > 5).all()
                if all_greater_than_5:
                    print(row['gene_id'],START,END,"Found in Father")
                    Final_DF = pd.concat([Final_DF, EXON_df], ignore_index=True)
                else:
                    print(row['gene_id'],START,END,"  Not Found in Father Either")
Final_DF.to_csv('Exons_Passed.tsv', sep='\t',index=False)
