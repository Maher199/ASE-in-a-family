import os,sys
import pandas as pd
from argparse import ArgumentParser

parser = ArgumentParser(description="This script is a part of the ASE pipeline to define the expression cases based on the log2fc and padj value after performing tests and comparison in family, it is designed for family from parents and  8 children" +
  "                                                                                                                ------------------------------------------------------------------- " +                          "  **** Genomics and Bioinformatics lab, M.Najjar, 2024 ****")

parser.add_argument('--input_file','-i',help="the dataset after running contrast tests in the family")
parser.add_argument('--output_dir','-o',help="work directory",default="./")
parser.add_argument('--exons_filtered','-e',help="input file of the filtered exons based on the coverage")
args = parser.parse_args()

Work_Dir = args.output_dir +"/"

## selected log2FC with Padj

padj_log2fc =  pd.read_csv(args.input_file,sep='\t')
padj_log2fc['father_log2fc'] = 0
padj_log2fc['father_padj'] = 0
padj_log2fc = padj_log2fc.round(3)
log2fc_columns = padj_log2fc.columns[padj_log2fc.columns.str.contains('log2fc')]

# Round the log2fc columns to 1 decimal
padj_log2fc[log2fc_columns] = padj_log2fc[log2fc_columns].round(1)

## The genes that passed the Exon filtering process, so we only keep the most relibale cases. This argument is optional
if args.exons_filtered:
    Genes_Filtered_by_Coverage = pd.read_csv(args.exons_filtered, sep='\t')

    ### Extract only the genes that passed the coverage

    padj_log2fc = padj_log2fc[padj_log2fc['Gene'].isin(Genes_Filtered_by_Coverage['Gene'])]

## H_L
aligned_children, aligned_mother = padj_log2fc[['child_1_log2fc', 'child_2_log2fc', 'child_3_log2fc', 'child_4_log2fc', 'child_5_log2fc', 'child_6_log2fc', 'child_7_log2fc', 'child_8_log2fc']].align(padj_log2fc[['.mother_log2fc']], axis=0, copy=False)

# Filtering DataFrames based on conditions
H_L = padj_log2fc[(padj_log2fc['.mother_log2fc'] > 1) & (padj_log2fc['.mother_padj'] < 0.05) & (aligned_mother.values >= aligned_children - 0.2).all(axis=1) & (aligned_children >= -0.2).all(axis=1)]



H_L.reset_index(drop=True, inplace=True)

H_L = H_L.drop_duplicates(subset=['Gene'])
H_L = H_L[
    ~(H_L[['child_1_log2fc', 'child_2_log2fc', 'child_3_log2fc', 'child_4_log2fc', 'child_5_log2fc', 'child_6_log2fc', 'child_7_log2fc', 'child_8_log2fc']].gt(1).any(axis=1) &
      H_L[['child_1_log2fc', 'child_2_log2fc', 'child_3_log2fc', 'child_4_log2fc', 'child_5_log2fc', 'child_6_log2fc', 'child_7_log2fc', 'child_8_log2fc']].le(0.4).any(axis=1))
]
H_L.to_csv(Work_Dir+"H_L.tsv",sep='\t',index=False)

### Create Symblic Dataset from the H_L log2FC, replace high expression with H, low expression with L and Moderate Expression with M

H_L_Symb = H_L.copy()
H_L_Symb = H_L_Symb[['Gene','.mother_log2fc','child_1_log2fc', 'child_2_log2fc', 'child_3_log2fc', 'child_4_log2fc', 'child_5_log2fc', 'child_6_log2fc', 'child_7_log2fc', 'child_8_log2fc','father_log2fc']]
H_L_Symb.columns = H_L_Symb.columns.str.replace('_log2fc$', '', regex=True)
H_L_Symb.rename(columns={'.mother': 'mother'}, inplace=True)

def categorize_value(val, column_name):
    if column_name == "mother":
        return 'H'
    elif column_name == "father":
            return 'L'
    else:
            return 'M'
for column in H_L_Symb.columns:
    if column != "Gene":
        H_L_Symb[column] = H_L_Symb[column].apply(lambda x: categorize_value(x, column))

H_L_Symb.to_csv(Work_Dir+"H_L_Symbolised.tsv",sep='\t',index=False)

## L_H
aligned_children, aligned_mother = padj_log2fc[['child_1_log2fc', 'child_2_log2fc', 'child_3_log2fc', 'child_4_log2fc', 'child_5_log2fc', 'child_6_log2fc', 'child_7_log2fc', 'child_8_log2fc']].align(padj_log2fc[['.mother_log2fc']], axis=0, copy=False)

# Filtering DataFrames based on conditions
L_H = padj_log2fc[(padj_log2fc['.mother_log2fc'] < -1) & (padj_log2fc['.mother_padj'] <= 0.05) & (aligned_mother.values <= aligned_children + 0.2).all(axis=1) & (aligned_children <= 0.2).all(axis=1)]



L_H.reset_index(drop=True, inplace=True)

L_H = L_H.drop_duplicates(subset=['Gene'])
L_H = L_H[
    ~(L_H[['child_1_log2fc', 'child_2_log2fc', 'child_3_log2fc', 'child_4_log2fc', 'child_5_log2fc', 'child_6_log2fc', 'child_7_log2fc', 'child_8_log2fc']].le(-1).any(axis=1) &
      L_H[['child_1_log2fc', 'child_2_log2fc', 'child_3_log2fc', 'child_4_log2fc', 'child_5_log2fc', 'child_6_log2fc', 'child_7_log2fc', 'child_8_log2fc']].gt(-0.4).any(axis=1))
]
L_H.to_csv(Work_Dir+"L_H.tsv",sep='\t',index=False)

## create Symbolic dataset for L_H

L_H_Symb = L_H.copy()
L_H_Symb = L_H_Symb[['Gene','.mother_log2fc','child_1_log2fc', 'child_2_log2fc', 'child_3_log2fc', 'child_4_log2fc', 'child_5_log2fc', 'child_6_log2fc', 'child_7_log2fc', 'child_8_log2fc','father_log2fc']]
L_H_Symb.columns = L_H_Symb.columns.str.replace('_log2fc$', '', regex=True)
L_H_Symb.rename(columns={'.mother': 'mother'}, inplace=True)

def categorize_value(val, column_name):
    if column_name == "mother":
        return 'L'
    elif column_name == "father":
            return 'H'
    else:
            return 'M'
for column in L_H_Symb.columns:
    if column != "Gene":
        L_H_Symb[column] = L_H_Symb[column].apply(lambda x: categorize_value(x, column))
L_H_Symb.to_csv(Work_Dir+"L_H_Symbolised.tsv",sep='\t',index=False)

## Define M_M cases
children_columns = ['child_1_log2fc', 'child_2_log2fc', 'child_3_log2fc', 'child_4_log2fc', 'child_5_log2fc', 'child_6_log2fc', 'child_7_log2fc', 'child_8_log2fc']
mother_threshold = 0.4
high_threshold = 0.8
# Filterinf for  cases where mother's expression is close to 0
close_to_zero_mother = padj_log2fc[padj_log2fc['.mother_log2fc'].abs() <= mother_threshold]

filtered_rows = []

for index, row in close_to_zero_mother.iterrows():
    # Check if at least one child's expression is higher than the threshold
    higher_than_mother = any(row[child_column] >= high_threshold for child_column in children_columns)
    # Check if at least one child's expression is lower than the mother's expression
    lower_than_mother = any(row[child_column] <= - high_threshold for child_column in children_columns)
    ## moderate Expression
    moderate = any(row[child_column] <= mother_threshold and row[child_column] >= - mother_threshold for child_column in children_columns)
    # If both conditions are met, append the row to the list
    if (higher_than_mother and lower_than_mother) or (higher_than_mother and moderate) or (lower_than_mother and moderate):
        filtered_rows.append(row)

M_M = pd.DataFrame(filtered_rows)
#df = close_to_zero_mother
range_min = -0.4
range_max = 0.4

mask1_lst = []
mask2_lst = []

for i in range(1, 9):  # Assuming there are 8 child columns (child_1 to child_8)
    log2fc_col = f'child_{i}_log2fc'
    padj_col = f'child_{i}_padj'
    
    # Check if log2fc > 1 and padj < 0.05
    mask1 = (M_M[log2fc_col] >= 1) & (M_M[padj_col] >= 0.05)
    mask2 = (M_M[log2fc_col] <= -1) & (M_M[padj_col] >= 0.05)
    mask1_lst.append(mask1)
    mask2_lst.append(mask2)
# Combine all masks
mask_1 = pd.concat(mask1_lst, axis=1).any(axis=1)
mask_2 = pd.concat(mask2_lst, axis=1).any(axis=1)
filtered_df_4 = M_M[mask_1]
filtered_df_5 = M_M[mask_2]
con = pd.concat([filtered_df_4, filtered_df_5])
con = con.drop_duplicates()
M_M = M_M[~M_M['Gene'].isin(con['Gene'])]

M_M.to_csv(Work_Dir+"M_M.tsv",sep='\t',index=False)

## Symbolise M_M

M_M_Symb = M_M.copy()
M_M_Symb = M_M_Symb[['Gene','.mother_log2fc','child_1_log2fc', 'child_2_log2fc', 'child_3_log2fc', 'child_4_log2fc', 'child_5_log2fc', 'child_6_log2fc', 'child_7_log2fc', 'child_8_log2fc','father_log2fc']]
M_M_Symb.columns = M_M_Symb.columns.str.replace('_log2fc$', '', regex=True)
M_M_Symb.rename(columns={'.mother': 'mother'}, inplace=True)

def categorize_value(val, column_name):
    if column_name == "mother":
        return 'M'
    elif column_name == "father":
            return 'M'
    else:
        if val >= 0.8:
            return 'H'
        elif val <= -0.8:
            return 'L'
        else:
            return 'M'
for column in M_M_Symb.columns:
    if column != "Gene":
        M_M_Symb[column] = M_M_Symb[column].apply(lambda x: categorize_value(x, column))

M_M_Symb.to_csv(Work_Dir+"M_M_Symbolised.tsv",sep='\t',index=False)

## H * M cases 

high_threshold = 0.8
low_threshold = 0.4
df = padj_log2fc.copy()
# Filtering for cases where the mother and some children are high in positive value
filtered_rows = []

for index, row in df.iterrows():
    if row['.mother_log2fc'] >= high_threshold:
        high_children = [child for child in children_columns if row[child] >= high_threshold]
        low_children = [child for child in children_columns if (row[child] >= - low_threshold) and ( row[child] <= low_threshold)]
        above_one = [child for child in children_columns if row[child] >= 1]
        if high_children and low_children and above_one:
            # Append the row to the filtered rows list
            filtered_rows.append(row)

H_M = pd.DataFrame(filtered_rows)

## filtering based on padj
mask1_lst = []
for i in range(1, 9):  # Assuming there are 8 child columns (child_1 to child_8)
    log2fc_col = f'child_{i}_log2fc'
    padj_col = f'child_{i}_padj'
    
    # Check if log2fc > 1 and padj < 0.5
    mask1 = (H_M[log2fc_col] >= 1) & (H_M[padj_col] >= 0.05)
    mask1_lst.append(mask1)
# Combine masks
mask_1 = pd.concat(mask1_lst, axis=1).any(axis=1)
filtered_ = H_M[mask_1]
H_M = H_M[~H_M['Gene'].isin(filtered_['Gene'])]



H_M.to_csv(Work_Dir+"H_M.tsv",sep='\t',index=False)

## Symbolise H_M

H_M_Symb = H_M.copy()
H_M_Symb = H_M_Symb[['Gene','.mother_log2fc','child_1_log2fc', 'child_2_log2fc', 'child_3_log2fc', 'child_4_log2fc', 'child_5_log2fc', 'child_6_log2fc', 'child_7_log2fc', 'child_8_log2fc','father_log2fc']]
H_M_Symb.columns = H_M_Symb.columns.str.replace('_log2fc$', '', regex=True)
H_M_Symb.rename(columns={'.mother': 'mother'}, inplace=True)

def categorize_value(val, column_name):
    if column_name == "mother":
        return 'H'
    elif column_name == "father":
            return 'M'
    else:
        if val >= 0.8:
            return 'H'

        else:
            return 'M'
for column in H_M_Symb.columns:
    if column != "Gene":
        H_M_Symb[column] = H_M_Symb[column].apply(lambda x: categorize_value(x, column))

H_M_Symb.to_csv(Work_Dir+"H_M_Symbolised.tsv",sep='\t',index=False)

## L * M cases 


high_threshold = 0.8
low_threshold = 0.4
df = padj_log2fc.copy()
### Filtering for cases where the mother and some children are high in positive value
filtered_rows = []

# Loop through the DataFrame row by row
for index, row in df.iterrows():
    if row['.mother_log2fc'] <= - high_threshold:
        high_children = [child for child in children_columns if row[child] <= - high_threshold]
        low_children = [child for child in children_columns if (row[child] >= - low_threshold) and ( row[child] <= low_threshold)]
        above_one = [child for child in children_columns if row[child] <= -1]
        if high_children and low_children and above_one:
            filtered_rows.append(row)

L_M = pd.DataFrame(filtered_rows)

mask1_lst = []
for i in range(1, 9):  # Assuming there are 8 child columns (child_1 to child_8)
    log2fc_col = f'child_{i}_log2fc'
    padj_col = f'child_{i}_padj'
    
    # Check if log2fc > 1 and padj < 0.05
    mask1 = (L_M[log2fc_col] <= -1) & (L_M[padj_col] >= 0.05)
    mask1_lst.append(mask1)
# Combine all masks to filter the DataFrame
mask_1 = pd.concat(mask1_lst, axis=1).any(axis=1)
filtered_ = L_M[mask_1]
L_M = L_M[~L_M['Gene'].isin(filtered_['Gene'])]


##Save
L_M.to_csv(Work_Dir+"L_M.tsv",sep='\t',index=False)

## Symbolise L_M

L_M_Symb = L_M.copy()
L_M_Symb = L_M_Symb[['Gene','.mother_log2fc','child_1_log2fc', 'child_2_log2fc', 'child_3_log2fc', 'child_4_log2fc', 'child_5_log2fc', 'child_6_log2fc', 'child_7_log2fc', 'child_8_log2fc','father_log2fc']]
L_M_Symb.columns = L_M_Symb.columns.str.replace('_log2fc$', '', regex=True)
L_M_Symb.rename(columns={'.mother': 'mother'}, inplace=True)

def categorize_value(val, column_name):
    if column_name == "mother":
        return 'L'
    elif column_name == "father":
            return 'M'
    else:

        if val <= -0.8:
            return 'L'
        else:
            return 'M'
for column in L_M_Symb.columns:
    if column != "Gene":
        L_M_Symb[column] = L_M_Symb[column].apply(lambda x: categorize_value(x, column))

L_M_Symb.to_csv(Work_Dir+"L_M_Symbolised.tsv",sep='\t',index=False)
