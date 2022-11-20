import pandas as pd
from functools import partial


def remap_multi(hgnc_list, valstring):
    return ';'.join([val for val in valstring.split(';') if val in hgnc_list])


# files to fix plus reference gene file
cnv1 = "processed/cnv_ref_fixed_1.tsv"
cnv2 = "processed/cnv_ref_fixed_3.csv"
var1 = "processed/variant_ref_21.tsv"
var2 = "processed/variant_ref_32.tsv"
var3 = "processed/variants.txt"
cnv1_df = pd.read_csv(cnv1, sep='\t')
cnv2_df = pd.read_csv(cnv2, sep='\t')
var1_df = pd.read_csv(var1, sep='\t')
var2_df = pd.read_csv(var2, sep='\t')
var3_df = pd.read_csv(var3, sep='\t')

# get list of current gene names
genes = "/Users/bkamphaus/code/pret/seed_data/raw/genes/hgnc_complete_set.txt"
genes_df = pd.read_csv(genes, sep='\t')
hugos = list(genes_df['symbol'])

# remap cnv files
cnv1_df['Genes'] = cnv1_df['Genes'].apply(partial(remap_multi, hugos))
cnv2_df['Genes'] = cnv2_df['Genes'].apply(partial(remap_multi, hugos))

# remove from var ref files genes not in hugos
keep_var1_df = var1_df[~((var1_df['variable'] == 'Hugo_Symbol') & (var1_df['value'].apply(lambda x: x not in hugos)))]
keep_var2_df = var2_df[~((var2_df['variable'] == 'Hugo_Symbol') & (var2_df['value'].apply(lambda x: x not in hugos)))]

# var3 in a different format
# drop if not here
# fix: var3_df['hugo'] = cnv2_df['hugo']

# write fixed files
keep_var1_df.to_csv(var1, sep='\t', header=True)
keep_var2_df.to_csv(var2, sep='\t', header=True)
cnv1_df.to_csv(cnv1, sep='\t', header=True)
cnv2_df.to_csv(cnv2, sep='\t', header=True)
