from scipy.stats import hypergeom  # hypergeometric distribution is a discrete
# distribution that models the number of events in a fixed sample size when you
# know the total number of items in the population that the sample is from
import pandas as p  # data manipulation and analysis
import json  # encode/decode JSON data
import math
import os  # functions for interacting with the operating system
import xlsxwriter

## Global paths
# json_dir contains a bunch of JSON files for GO CC, BP, and MF (combined and separated)
json_dir = 'C:\\Users\\peree\\OneDrive\\Desktop\\CompBio_Code\\fun_enrich_files\\mapping_files'
# genefile contains all yeast ORFS (6626) mapped to gene name OR themselves if no gene name is available
genefile = os.path.join(json_dir, 'gene_map_20210527.json')

# read file and load in genes
with open(genefile, 'r') as fp:
    genes = json.load(fp)

# Global variables
cluster_dir = 'C:\\Users\\peree\\OneDrive\\Desktop\\CompBio_Code\\fun_enrich_files'
out_dir = 'C:\\Users\\peree\\OneDrive\\Desktop\\CompBio_Code\\fun_enrich_files\\outputs'

split_CC = ['no_cell_cycle', 'with_cell_cycle']
split_AS = ['with_AreaShape', 'without_AreaShape']
split_genes = ['all_genes', 'essential_only', 'no_vesicle_genes', 'non_essential_only']

# Probably need to change these depending on the datasets under study
markers = ['Cdc11', 'Dad2', 'Heh2', 'Hta2', 'Nop10', 'Nuf2', 'Om45', 'Pil1',
           'Psr1', 'Rad52', 'Sec21', 'Sec7', 'Spf1']

# JSON files
go_slim_ontologies = ['go_slim_c_20210430', 'go_slim_p_20210430']

# Dictionary for remapping column names?
col_map = {'go_id': 'GO ID', 'go_name': 'Ontology', 'pval': 'P-value', 'hits_in_term': 'Hits in term',
           'term_size': 'Term size', 'all_hits': 'All hits', 'all_background': 'All background',
           'genes_orf': 'Genes (ORF)', 'genes_name': 'Genes (Standard Name)',
           'bonferroni': 'P-value (Bonferroni)', 'fold_enrichment': 'Fold enrichment'}
col_final = ['GO ID', 'Ontology', 'P-value (Bonferroni)', 'P-value', 'Hits in term', 'Term size', 'All hits',
            'All background', 'Fold enrichment', 'Genes (ORF)', 'Genes (Standard Name)']


# from thecellmap (Matej Usaj)
# annot is a dictionary of GO terms (keys) and genes mapped to these terms?
def do_enrichment_analysis(query, background, annot):
    # Background is all yeast ORFs and query is list of genes in i.e., each cluster
    m = len(background)  # M
    nn = len(query)  # N

    vals = []
    max_hits = 0

    # For each term (key) and genes mapped to it...
    for term, term_genes in annot.items():
        # background is a Set; category contains everything that's common between
        # background and term_genes for a given annotation
        category = background.intersection(term_genes)
        n = len(category)  # how many genes in background belong to a particular annotation
        hits = query.intersection(category)  # genes common between query list and category
        x = len(hits)  # how many hits?
        max_hits = max(max_hits, x)  # once loop finishes, get the max number of hits across all annotations

        if 0 in (n, x):
            # either no genes in this term or no hits in this term
            continue

        # use hypergeometric distribution to calculate p-value
        vals.append(
            term + (hypergeom.sf(x - 1, m, n, nn), x, n, nn, m, ','.join(hits)))

    # Start creating the table that will be part of the final output
    df = p.DataFrame(vals,
                     columns=['go_id', 'go_name', 'pval', 'hits_in_term',
                              'term_size', 'all_hits',
                              'all_background', 'genes_orf'])
    df = df.sort_values('pval')  # sort p-values from least to greatest
    df.loc[:, 'bonferroni'] = df.pval * df.shape[0]  # add new bonferroni column, calculate it as p-val * nrows?
    df.loc[:, 'bonferroni'] = df.loc[:, 'bonferroni'].clip(upper=1)  # force any number that's >1 to be changed to 1 (numbers <1 are unchanged)
    df.loc[:, 'fold_enrichment'] = (df.hits_in_term / df.all_hits) / \
                                   (df.term_size / df.all_background)  # add new fold enrichment column

    df = df.loc[df.hits_in_term > 0]  # filter out rows where hits_in_term is 0

    # add gene names
    genes_name = []
    # for every ORF in genes_orf column...
    for orfs in df['genes_orf'].values:
        orf_to_name = []
        for orf in orfs.split(','):  # split orfs separated by a comma into a list, then iterate on each orf
            orf_to_name.append(genes[orf])  # genes is a dictionary with orfs as keys and gene names assigned to these keys
                                            # append the gene name by searching genes dictionary using orf key
        genes_name.append(','.join(orf_to_name))  # add contents of orf_to_name list to genes_name list
    df['genes_name'] = genes_name  # add new column
    return df

# load annotation files
all_slims = {}
for ontology in go_slim_ontologies:  # for JSON file in list of JSON files...
    ontology_jsonfile = os.path.join(json_dir, '%s.json' % ontology)  # output looks weird (has \\)

    with open(ontology_jsonfile, 'r') as fp:
        annotation_dict = json.load(fp)  # converts JSON file to iterable python object
        # the object consists of nested dictionaries and lists

    annotation = list(annotation_dict.keys())[0]  # get the topmost key from this JSON object
    slim = {}
    terms = []
    for term, values in annotation_dict[annotation].items():
        # make list of all possible terms in the GO Slim file
        # add new entry to slim dict: {(GO_ID, term_name): [ORFs]}
        if term not in terms:
            terms.append(term)
            slim[tuple(values[0])] = values[1]
    all_slims[ontology] = slim  # add contents of slim dictionary to all_slims under new key
    print("%s — number of terms: %d" % (annotation, len(terms)))

all_slims.keys()

for CC in split_CC:
    for AS in split_AS:
        for genes_dir in split_genes:
            # load in input, query, and background
            ds_dir = '%s\\%s\\%s' % (CC, AS, genes_dir)  # marker_rep_PCAdataset
            strain_path = os.path.join(cluster_dir, '%s\\input.csv' % ds_dir)
            strain_df = p.read_csv(strain_path)

            query_path = os.path.join(cluster_dir, '%s\\query.csv' % ds_dir)
            query = p.read_csv(query_path)
            query = set(query['ORF'])

            background_path = 'C:\\Users\\peree\\OneDrive\\Desktop\\CompBio_Code\\fun_enrich_files\\genes_list\\deletion_and_ts_array_genes.csv'
            background = p.read_csv(background_path)
            background = set(background['ORF'])
            len_background = len(background)
            N = round(len_background*0.01)  # 1% of number of unique ORFs in background (??? I still don't get it)
            print(f'{CC}/{AS}/{genes_dir} ‒ {strain_df.shape} strain_df: {strain_df.shape[0]} ‒ {len_background} unique ORFs ‒ N: {N}')

            # get cluster column names for all cell cycle stages
            stages = []
            columns = strain_df.columns
            for c in columns:
                if c not in ['ORF', 'Name', 'Allele', 'Strain ID']:
                    stages.append(c)

            # go through each annotation
            for ontology in go_slim_ontologies:
                slim = all_slims[ontology]
                ontology_name = '_'.join(ontology.split('_')[:-1])  # 'mc_slim_20180507' -> ['mc', 'slim', '20180507'] -> 'mc_slim'

                out_path = os.path.join(out_dir, '%s-%s-%s-%s.xlsx' % (CC, AS, genes_dir, ontology_name))
                writer = p.ExcelWriter(out_path, engine='xlsxwriter')
                print(os.path.relpath(out_path, start=cluster_dir))  # get relative output path from cluster_dir

                # go through each appended cluster, make a copy of strain_df, and sort from greatest to least
                for stage in stages:
                    sheetname = stage.replace('/', '-')
                    try:
                        df_stage = strain_df.sort_values(by=stage, ascending=False)
                    except KeyError:
                        continue

                    # do enrichment analysis
                    df_slim = do_enrichment_analysis(query, background, slim)
                    df_slim.rename(columns=col_map, inplace=True)
                    df_slim = df_slim[col_final]
                    df_slim.to_excel(writer, sheet_name=sheetname, index=False)
                writer.close()

                # query = set()
                # for sid, straindata in df_stage.iterrows():
                #    orf = straindata['ORF']
                #    stage_val = straindata[stage]
                #    if orf not in query and stage_val > 0.05:
                #        query.add(orf)
print("DONE")
