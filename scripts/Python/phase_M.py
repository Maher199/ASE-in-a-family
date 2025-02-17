#~/miniconda3/bin/python3


# ------- Genomics and Bioinformatics lab,  M.Najjar, 2024 ------- #

# Phase_M phases the haplotype of the parent of a family whenever there are supporting variants in this region, heterozygous in one parent and homozygous in the other parent

## It checks for the abundancy of the variants pattern in the given region 

## The programme needs pandas and matpltlib libraries to be installed

import pandas as pd
import os,sys
import matplotlib.pyplot as plt
from matplotlib import rcParams
from argparse import ArgumentParser

parser = ArgumentParser(description="Phase_M phases the haplotype of the parent of a family whenever there are supporting variants in this region, heterozygous in one parent and homozygous in the other parent" +
  "                                                                                                                ------------------------------------------------------------------- " +                                                                            "                                     ------------------------------------------------------------------- " +                          "  **** Genomics and Bioinformatics lab,  M.Najjar, 2024 ****")

parser.add_argument('--input_file','-i',help="The input VCF file")
parser.add_argument('--output_file','-o',help="The Final phased image",default="Haplotype.png")
parser.add_argument('--chromosme','-chr',help="Specify the chromosme needed to be phased")
parser.add_argument('--start','-s',type=int,help="Start region to be phased",)
parser.add_argument('--end','-e',type=int,help="End region to be phased")
args = parser.parse_args()


### Parsing VCF and extract the corresponding region

def parse_vcf(vcf,Chr,Start,End):
    vcf_file = pd.read_csv(vcf, comment="#", header=None, sep="\t")
    columns_names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO','INFO_GT','mother', 'child_1', 'child_2', 'child_3', 'child_4', 'child_5', 'child_6', 'child_7', 'child_8', 'father']


    vcf_file.columns = columns_names
    genotype_df = pd.DataFrame()

    # Add chromosome and position columns
    genotype_df["CHROM"] = vcf_file["CHROM"]
    genotype_df["POS"] = vcf_file["POS"]

    individuals = vcf_file.columns[9:]  # Assuming individual columns start from index 9
    for ind in individuals:
        genotype_df[f"{ind}"] = vcf_file[ind].apply(lambda x: x.split(":")[0])
    Extracted = genotype_df[(genotype_df['CHROM'] == Chr) & (genotype_df['POS'] >= Start) & (genotype_df['POS'] <= End)]
    Father_HET = Extracted[(Extracted['mother'].isin(['1|1','1/1','0|0','0/0'])) & (Extracted['father'].isin(['0/1','0|1']))]
    Mother_HET = Extracted[(Extracted['father'].isin(['1|1','1/1','0|0','0/0'])) & (Extracted['mother'].isin(['0/1','0|1']))]

    return Father_HET, Mother_HET

vcf_file=args.input_file
filename = args.output_file
Chr = args.chromosme
Start = args.start
End = args.end

Father_HET, Mother_HET = parse_vcf(vcf_file,Chr,Start,End)

### define the grouping of the children (like the parent, others)

children = set()

def find_grouping(DF,parent):
    GROUOS = {}
    global children
    for _,ROW in DF.iterrows():
        Parent_like = []
        for IND in DF.columns:
            if IND.startswith("child"):
                children.add(IND)
                if ROW[IND] == ROW[parent]:
                    Parent_like.append(IND)
        Parent_like = tuple(Parent_like)
        if Parent_like in GROUOS:
            GROUOS[Parent_like] += 1
        else:
            GROUOS[Parent_like] = 1
    Abundance = max(GROUOS,key=GROUOS.get)
    return Abundance
mother_like_hap = find_grouping(Mother_HET,'mother')
father_like_hap = find_grouping(Father_HET,'father')

### combine the parents' haplotypes and sort the DF for plotting

column_names = ['sample', 'seqname', 'start', 'end', 'haplotype']


FINAL_DF = pd.DataFrame(columns=column_names)

#define a function to append to the Dataframe
def append_to_final_df(child, hap_type, hap_label):
    haplotype = f"{hap_type}{hap_label}"
    row = [child, Chr, Start, End, haplotype]
    return pd.DataFrame([row], columns=column_names)

# Loop over children and assign parents' haplotypes
for CHILD in children:
    if CHILD in mother_like_hap:
        FINAL_DF = pd.concat([FINAL_DF, append_to_final_df(f"{CHILD}_M", "M", 1)], ignore_index=True)
    else:
        FINAL_DF = pd.concat([FINAL_DF, append_to_final_df(f"{CHILD}_M", "M", 2)], ignore_index=True)
    
    if CHILD in father_like_hap:
        FINAL_DF = pd.concat([FINAL_DF, append_to_final_df(f"{CHILD}_F", "F", 1)], ignore_index=True)
    else:
        FINAL_DF = pd.concat([FINAL_DF, append_to_final_df(f"{CHILD}_F", "F", 2)], ignore_index=True)



###PLOTTING

### parameters here should be changed to accomodiate the region of interest,e.g. colors and hashes

import matplotlib.patches as patches
#define haplotype colors and hatch 
haplotype_colors = {'M1': 'red', 'M2': 'green', 'F1': 'red', 'F2': 'green'}
haplotype_hatches = {'M1': '', 'M2': '', 'F1': '///', 'F2': '///'}
#haplotype_labels = {'M1': '    B', 'M2': '    B', 'F1': '    A', 'F2': '    B'}  # Letters to plot on bars
df = FINAL_DF.copy()
def plot_haplotype_distribution(df, region_start, region_end, filename):
    plt.figure(figsize=(17, 13))
    
    desired_order = ['child_1_M', 'child_1_F', 'child_4_M', 'child_4_F', 
                 'child_7_M', 'child_7_F', 'child_3_M', 'child_3_F',
                 'child_5_M', 'child_5_F', 'child_8_M', 'child_8_F',
                 'child_2_M', 'child_2_F', 'child_6_M', 'child_6_F']

    # Convert 'sample' column to categorical
    df['sample'] = pd.Categorical(df['sample'], categories=desired_order, ordered=True)

    # Sort the Dataframe by 'sample' 
    df_sorted = df.sort_values('sample').reset_index(drop=True)

    
    df_sorted['order'] = df.index
    df = df_sorted.copy()
    # Define gap sizes
    gap_within_pairs = 1
    gap_between_pairs = 2.0
    
    # Positioning
    num_samples = len(df)
    x_positions = []
    
    current_position = 0
    for idx in range(num_samples):
        x_positions.append(current_position)
        if idx % 2 == 1:  # After every pair
            current_position += gap_between_pairs
        else:  
            current_position += gap_within_pairs
    
    for idx, (x_pos, row) in enumerate(zip(x_positions, df.itertuples(index=False))):
        if row.end > region_start and row.start < region_end:
            bar = plt.bar(
                x_pos, 
                min(row.end, region_end) - max(row.start, region_start), 
                bottom=max(row.start, region_start), 
                color=haplotype_colors[row.haplotype], 
                hatch=haplotype_hatches[row.haplotype],
                edgecolor='black',  # Show haatches
                label=row.haplotype,
                width=0.8 
            )
            plt.rcParams['hatch.linewidth'] = 2.0
            # Position the text at the beginnin
            text_y = max(row.start, region_start) + 5  # Adding a small offset
            bar_center_x = bar[0].get_x() + bar[0].get_width() / 2

            ## Add text if nesseccary
            #plt.text(bar_center_x, text_y, haplotype_labels[row.haplotype], 
                     #ha='center', va='bottom', fontsize=15, color='white', weight='extra bold', rotation=90)
    plt.ylim(region_start, region_end)
    plt.ylabel('Position')
    plt.xticks(x_positions, df['sample'], rotation=90) 
    plt.xlabel('Sample')
    plt.title('Haplotype Distribution')
    
    # Add legend, place it 
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(sorted(zip(labels, handles), key=lambda x: x[0]))
    plt.legend(by_label.values(), by_label.keys(), title="Haplotype", loc='lower right')

    # Save
    plt.savefig(filename, dpi=300, bbox_inches='tight')  #high resolution

    plt.show()

# Plotting
plot_haplotype_distribution(df, Start, End, filename)
