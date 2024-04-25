"""
MAIN GOAL: Create a selectio for the genes to test + plots:
-- Intersection gene in at least one GEL gene-panel and essential in HAP1 -> Initial number of genes: 865
-- genes with higher GEL status = 3 are kepte
-- Number of genes in current list: 136
The 136 genes and summary information are safe as well as the initial 865 genes with different filters as boolean cols.
--> this is the old list now all genes are used and they are selected later on in variant_selection_control
To update from ENSG 37 to 38
ZNHIT3 (zinc finger HIT-type containing 3)
EnsemblGeneIds (GRCh38): ENSG00000273611
EnsemblGeneIds (GRCh37): ENSG00000108278

RBM8A (RNA binding motif protein 8A)
EnsemblGeneIds (GRCh38): ENSG00000265241
EnsemblGeneIds (GRCh37): ENSG00000131795

Date: 14/02/2024
Author: Chloé Terwagne
"""
# IMPORT ---------------------------------------------------------------------------------------------------------------

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
import datetime
import urllib.request
import plotly.express as px
from collections import Counter
import json
import sys

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 2000)

# CONSTANT ------------------------------------------------------------------------------------------------------------
DATE = str(datetime.date.today()).replace('-', '')
PLOT = False  # To indicate whether to plot or not
SAVING_PLOT = False  # To indicate whether to save the plot or show
SAVING_DF = False  # To indicate whether to save the dataframe or not

# Path to save figures
FIG_PATH = 'data/figures/' + str(datetime.date.today()).replace('-', '') + '_'

# DF ----------------------------------------------------------------------------------------------------------------
# Blomen et al., 2015
HAP1_blomen_df = pd.read_csv("/Users/terwagc/PycharmProjects/EG_visualisation/data/input_df/aac7557_sm_table_s2_blomen_2015.csv", header=1)

# core essentiality in human cell line in OGEEv3 - https://v3.ogee.info/#/downloads
df_ogee = pd.read_csv('/Users/terwagc/PycharmProjects/EG_visualisation/data/input_df/CSEGs_CEGs.txt', sep='\t', skiprows=11)

# Gnomad O/E -- v4
df_oe_all = pd.read_csv('/Users/terwagc/PycharmProjects/EG_visualisation/data/input_df/release_v4.0_constraint_gnomad.v4.0.constraint_metrics.tsv',
                        sep='\t')
# gene panel
panel_df = pd.read_csv('/Users/terwagc/PycharmProjects/EG_visualisation/data/input_df/panel_app_out/bigbedtool/genesv2.bed', sep='\t', header=None)

# Position
gene_gtf_df = pd.read_csv('/Users/terwagc/PycharmProjects/EG_visualisation/data/input_df/gtf_genes_position_files/Homo_sapiens.GRCh38.110_genes.txt',
                          sep='\t')

# ClinVar gene summary
file_url = 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/gene_specific_summary.txt'
destination_path = '/Users/terwagc/PycharmProjects/EG_visualisation/data/input_df/' \
                   + str(datetime.date.today()).replace('-', '') + '_ClinVar_gene_specific_summary.txt'
urllib.request.urlretrieve(file_url, destination_path)
clinvar_gene_specific_summary = pd.read_csv(destination_path, sep='\t', header=1)


# FUNCTION ------------------------------------------------------------------------------------------------------------


def cleaning_panel_bed(df):
    columns = ['chrom_hg19', 'start_hg19', "end_hg19", 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb',
               'blockCount', 'blockSizes',
               'blockStarts', 'geneSymbol', 'biotype', 'hgncID', 'geneName', 'omimGene', 'ensemblGenes', 'entityType',
               'entityName',
               'confidenceLevel', 'penetrance', 'modeOfPanthogenicity', 'publications', 'evidence', 'phenotypes',
               'modeOfInheritance',
               'tags', 'panelID', 'panelName', 'diseaseGroup', 'diseaseSubgroup', 'status', 'panelVersion',
               'versionCreated',
               'relevantDisorders', 'mouseOverField']

    # Assign columns to the DataFrame
    df.columns = columns

    # one row got split into two due to a missing comma between two publication ID-> Merging row 21534 and 21535 together
    pub1 = df.iloc[21534, 23]
    pub2 = df.iloc[21535, 0]
    df.iloc[21534, 23] = str(pub1) + ', ' + str(pub2)
    df.iloc[21534, 24:] = df.iloc[21535, 1:14].values
    df = df.drop(21535).reset_index(drop=True)
    df["chrom_hg19"] = df["chrom_hg19"].str.strip('chr')

    # Remove genes on the MT chromosome
    df = df[df['chrom_hg19'] != 'MT']

    # removing 37  coordinates
    df_panels = df[
        ['name', 'score', 'strand', 'blockCount', 'blockSizes', 'blockStarts', 'geneSymbol', 'biotype', 'hgncID',
         'geneName', 'omimGene', 'ensemblGenes', 'entityName', \
         'confidenceLevel', 'penetrance', 'modeOfPanthogenicity', 'publications', 'evidence', 'phenotypes',
         'modeOfInheritance', \
         'tags', 'panelID', 'panelName', 'diseaseGroup', 'diseaseSubgroup', 'status', 'panelVersion',
         'versionCreated', 'relevantDisorders', 'mouseOverField']]
    return df_panels


def dict_modeofinheritance_modeofpatho(df):
    # Create an empty dictionary to store the mode of inheritance and associated gene ENSEMBL_IDs
    mode_of_inheritance_dict = {}
    mode_of_pathogenecity_dict = {}

    # Define the mapping of current values to the desired shorter names
    mapping = {
        'BIALLELIC, autosomal or pseudoautosomal': 'BIALLELIC',
        'MONOALLELIC, autosomal or pseudoautosomal, imprinted status unknown': 'MONOALLELIC, unknown imprinting',
        'MONOALLELIC, autosomal or pseudoautosomal, NOT imprinted': 'MONOALLELIC, not imprinted',
        'BOTH monoallelic and biallelic, autosomal or pseudoautosomal': 'BOTH mono- and biallelic',
        'Unknown': 'Unknown',
        'BOTH monoallelic and biallelic (but BIALLELIC mutations cause a more SEVERE disease form), autosomal or pseudoautosomal': 'BOTH mono- and biallelic, severe',
        'Other': 'Other'
    }

    # Replace the values in the "modeOfInheritance" column
    df['modeOfInheritance'] = df['modeOfInheritance'].replace(mapping)
    # Iterate through the DataFrame rows
    for index, row in df.iterrows():
        # Extract the mode of inheritance and ENSEMBL_ID values from the row
        mode_of_patho = row['modeOfPanthogenicity']
        mode_of_inheritance = row['modeOfInheritance']
        ensembl_id = row['ENSEMBL_ID']

        # Update mode_of_inheritance_dict
        if ensembl_id in mode_of_inheritance_dict:
            mode_of_inheritance_dict[ensembl_id][mode_of_inheritance] += 1
        else:
            mode_of_inheritance_dict[ensembl_id] = Counter({mode_of_inheritance: 1})

        # Update mode_of_pathogenecity_dict
        if ensembl_id in mode_of_pathogenecity_dict:
            mode_of_pathogenecity_dict[ensembl_id][mode_of_patho] += 1
        else:
            mode_of_pathogenecity_dict[ensembl_id] = Counter({mode_of_patho: 1})

    return mode_of_inheritance_dict, mode_of_pathogenecity_dict


def cleaning_oe(df_oe_all):
    df_oe_all = df_oe_all[['gene', 'transcript', 'mane_select', 'lof.oe', 'mis.oe', 'syn.oe', 'lof.oe_ci.lower',
                           'mis.oe_ci.lower', 'syn.oe_ci.lower', 'lof.oe_ci.upper', 'mis.oe_ci.upper',
                           'syn.oe_ci.upper']]
    df_oe_all.rename(columns={'gene': 'GENE_SYMBOL'}, inplace=True)
    # Correcting annotation for POLR2A as  135423 POLR2A  NM_000937.5 is set to 'mane_select' == False
    df_oe_all.loc[135423, 'mane_select'] = True
    df_oe_all = df_oe_all[df_oe_all['mane_select'] == True]  # Get MANE transcript

    # All genes are duplicated with != transcript col (Ensembl and oRefSeq transcript ID but value are similar)
    df_oe_all = df_oe_all.sort_values(by=['GENE_SYMBOL', 'transcript'])
    df_oe_all.drop_duplicates("GENE_SYMBOL", inplace=True)  # keep only Ensembl transcript name
    return df_oe_all


def print_info_gene_panel(list_genes, df_gene_panel, gene_info):
    df_gene_panel_subset = df_gene_panel[df_gene_panel['ENSEMBL_ID'].isin(list_genes)]
    print('\nINFO PANELS FOR ' + str(len(list_genes)) + ' GENES')
    print(' Number of gene-panel pairs:', df_gene_panel_subset.shape[0])
    print(' Number of unique panels:', len(df_gene_panel_subset['panelID'].unique()))
    print(' Number of unique genes:', len(df_gene_panel_subset['geneSymbol_panel_app'].unique()))
    print(' Number gene-panel not Protein Coding:',
          len(df_gene_panel_subset[df_gene_panel_subset['biotype'] != 'Protein Coding']))

    # sort by GEL_Status to keep the highest Gel status
    df = df_gene_panel_subset.sort_values(['ENSEMBL_ID', "confidenceLevel"])
    df = df.drop_duplicates('ENSEMBL_ID', keep='last')
    green_df, orange_df, red_df = df[df.confidenceLevel == 3], df[df.confidenceLevel == 2], df[
        (df.confidenceLevel == 0) | (df.confidenceLevel == 1)]
    print(' Number "highest GEL status" for each unique gene across panels (red, amber, green)', red_df.shape[0],
          orange_df.shape[0], green_df.shape[0])

    for gene in list_genes:
        df_gene_panel_gene = df_gene_panel_subset[df_gene_panel_subset['geneSymbol_panel_app'] == gene]
        if gene_info:
            print('\tINFO FOR', gene)
            print('\t Number of unique panel associated to the gene:', len(df_gene_panel_gene['panelID'].unique()))
            print('\t Disease group:', list(df_gene_panel_gene['diseaseGroup']))
            print('\t panelName group:', list(df_gene_panel_gene['panelName']))
            print('\t confidenceLevel group:', list(df_gene_panel_gene['confidenceLevel']))
    return df_gene_panel_subset

def get_gene_higher_gel_status(gene_list, gene_panel_df):
    # Initialize an empty dictionary
    dict_gene_higher_confidence_level = {}
    for gene in gene_list:
        # Filter the DataFrame for the current gene
        gene_subset = gene_panel_df[gene_panel_df['ENSEMBL_ID'] == gene]
        # Find the maximum confidence level for the current gene
        max_confidence_level = gene_subset['confidenceLevel'].max()
        # Update the dictionary with the gene and its maximum confidence level
        dict_gene_higher_confidence_level[gene] = max_confidence_level
    return dict_gene_higher_confidence_level


def get_ensembl_list(ogee_df):
    def get_list(df, col):
        df = df.dropna(subset=[col])
        # Split the 'ensembl' column by ';'and flatten the nested list
        ensembl_lists = df[col].str.split(';')
        ensembl_flat_list = [item.strip() for sublist in ensembl_lists for item in sublist]
        return ensembl_flat_list

    df_ogee_CSEGs = ogee_df[ogee_df['essentiality'] == 'CSEGs']
    df_ogee_CEGs = ogee_df[ogee_df['essentiality'] == 'CEGs']

    ogee_CSEGs_ens_list = get_list(df_ogee_CSEGs, 'ensembl')
    ogee_CEGs_ens_list = get_list(df_ogee_CEGs, 'ensembl')

    ogee_CSEGs_gene_list = get_list(df_ogee_CSEGs, 'gene')
    ogee_CEGs_gene_list = get_list(df_ogee_CEGs, 'gene')
    # Note : there is NAN for ENSG
    return ogee_CSEGs_ens_list, ogee_CEGs_ens_list, ogee_CSEGs_gene_list, ogee_CEGs_gene_list


def get_list_subset_df(ratio_thr, clinvar_thr, df):
    df = df[(df['ratio'] <= ratio_thr) & (df['Alleles_reported_Pathogenic_Likely_pathogenic'] > clinvar_thr)]
    return df['ENSEMBL_ID'].to_list()


# PLOTS FUNCTION-----------------------------------------------------------------------------------------------------
def plt_swarm_essratio_highergelstatus(df_genes):
    # Create swarm plot
    plt.figure(figsize=(8, 6))  # Adjust the figure size if needed
    sns.swarmplot(data=df_genes, x='higher_gel_status', y='ratio',
                  palette={'1.0': 'red', '2.0': 'orange', '3.0': 'green'})

    # Add labels and title
    plt.xlabel('Higher GEL Status')
    plt.ylabel('Ratio')
    plt.title(
        'Gene ratio by higher GEL status across panel (N=' + str(df_genes.shape[0]) + ')')
    if SAVING_PLOT:
        plt.savefig(FIG_PATH + 'plt_swarm_essratio_highergelstatus.png')
    else:
        plt.show()


def plt_scatter_oe_bhap1ratio(df):
    # Define colors based on higher_gel_status
    color_map = {1: 'red', 2: 'orange', 3: 'green'}

    # Create subplots
    fig, axs = plt.subplots(3, 1, figsize=(20, 8), sharex=True)
    df_nonan = df.dropna(subset=['ratio', 'lof.oe', 'mis.oe', 'syn.oe'])
    print("number gene without NaN value: ", df_nonan.shape[0])
    # Iterate over each subplot
    for i, (ax, y_col) in enumerate(zip(axs, ['lof.oe', 'mis.oe', 'syn.oe'])):
        # Scatter plot for each subgroup
        for status, group in df_nonan.groupby('higher_gel_status'):
            ax.scatter(group['ratio'], group[y_col], color=color_map[status], label=f'Higher Gel Status {status}',
                       marker='o', alpha=0.5)

        # Calculate correlation coefficient and p-value
        correlation_coef, p_value = pearsonr(df_nonan['ratio'], df_nonan[y_col])

        ax.set_ylabel(y_col)
        ax.set_title(f'\nCorrelation Coef: {correlation_coef:.3f} P-value: {p_value:.5f}')
        ax.grid(True)
        ax.legend()

    # Set common X label
    axs[-1].set_xlabel('HAP1 Essentiality Ratio')

    # Title for the entire figure
    plt.suptitle('Observed over Expected mutation rate vs Essentiality ratio in HAP1 by variant type (N genes = ' + str(
        df_nonan.shape[0]) + ')',
                 y=0.98)

    if SAVING_PLOT:
        plt.savefig(FIG_PATH + 'plt_scatter_oe_bhap1ratio.png')
    else:
        plt.show()


def plt_hist_oe_by_ogeev3_ess(df_genes):
    # Filter rows with CSEGs == True, CEGs == True, and all rows
    df_csegs_true = df_genes[df_genes['CSEGs']]
    df_cegs_true = df_genes[df_genes['CEGs']]
    df_all = df_genes

    # Create the subplots with shared x-axis
    fig, axs = plt.subplots(3, 1, figsize=(18, 8), sharex=True)

    # Plot for lof.oe
    axs[0].hist(df_all['lof.oe'].dropna(), bins=30, alpha=0.1,
                label='all genes, N=' + str(df_all['lof.oe'].dropna().shape[0]))
    axs[0].hist(df_cegs_true['lof.oe'].dropna(), bins=30, alpha=0.5,
                label='Core-essential genes, N=' + str(df_cegs_true['lof.oe'].dropna().shape[0]))
    axs[0].hist(df_csegs_true['lof.oe'].dropna(), bins=30, alpha=0.8,
                label='Cancer-specific essential genes, N=' + str(df_csegs_true['lof.oe'].dropna().shape[0]))
    axs[0].set_title('Distribution of lof.oe Values')
    axs[0].set_ylabel('Number of genes')
    axs[0].legend()

    # Plot for mis.oe
    axs[1].hist(df_all['mis.oe'].dropna(), bins=30, alpha=0.1,
                label='all genes, N=' + str(df_all['mis.oe'].dropna().shape[0]))
    axs[1].hist(df_cegs_true['mis.oe'].dropna(), bins=30, alpha=0.5,
                label='Core-essential genes, N=' + str(df_cegs_true['mis.oe'].dropna().shape[0]))
    axs[1].hist(df_csegs_true['mis.oe'].dropna(), bins=30, alpha=0.8,
                label='Cancer-specific essential genes, N=' + str(df_csegs_true['mis.oe'].dropna().shape[0]))
    axs[1].set_title('Distribution of mis.oe Values')
    axs[1].set_ylabel('Number of genes')
    axs[1].legend()

    # Plot for syn.oe
    axs[2].hist(df_all['syn.oe'].dropna(), bins=30, alpha=0.1,
                label='all genes, N=' + str(df_all['syn.oe'].dropna().shape[0]))
    axs[2].hist(df_cegs_true['syn.oe'].dropna(), bins=30, alpha=0.5,
                label='Core-essential genes, N=' + str(df_cegs_true['syn.oe'].dropna().shape[0]))
    axs[2].hist(df_csegs_true['syn.oe'].dropna(), bins=30, alpha=0.8,
                label='Cancer-specific essential genes, N=' + str(df_csegs_true['syn.oe'].dropna().shape[0]))
    axs[2].set_title('Distribution of syn.oe Values')
    axs[2].set_xlabel('observed/expected ratio')
    axs[2].set_ylabel('Number of genes')
    axs[2].legend()

    # Title for the entire figure
    plt.suptitle('Histogram of the observed over expected mutation rate by variant type (N genes = ' + str(
        df_all['lof.oe'].dropna().shape[0]) + ')',
                 y=0.98)
    plt.tight_layout()

    if SAVING_PLOT:
        plt.savefig(FIG_PATH + 'plt_hist_oe_by_ogeev3_ess.png')
    else:
        plt.show()


def plt_hist_ratio_by_ogeev3_ess(df_genes):
    # Filter rows with CSEGs == True, CEGs == True, and all rows
    df_csegs_true = df_genes[df_genes['CSEGs']]
    df_cegs_true = df_genes[df_genes['CEGs']]
    df_all = df_genes

    # Create the plot
    plt.figure(figsize=(12, 6))

    # Plot
    plt.hist(df_all['ratio'].dropna(), bins=30, alpha=0.2,
             label='all genes, N=' + str(df_all['ratio'].dropna().shape[0]))
    plt.hist(df_cegs_true['ratio'].dropna(), bins=30, alpha=0.5,
             label='Core-essential genes, N=' + str(df_cegs_true['ratio'].dropna().shape[0]))
    plt.hist(df_csegs_true['ratio'].dropna(), bins=30, alpha=0.8,
             label='Cancer-specific essential genes, N=' + str(df_csegs_true['ratio'].dropna().shape[0]))

    # Add labels and title
    plt.xlabel('HAP1 Essentially Ratio')
    plt.ylabel('Number of genes')
    plt.title('Distribution of HAP1 essential ratio values by OGEE V3 essentially classification')
    plt.legend()

    if SAVING_PLOT:
        plt.savefig(FIG_PATH + 'plt_hist_ratio_by_ogeev3_ess.png')
    else:
        plt.show()


def plt_hist_clinvar_path(df, capped=60):
    # Create the plot
    plt.figure(figsize=(12, 6))

    # Replace '-' with NaN and convert to integers
    df['Alleles_reported_Pathogenic_Likely_pathogenic'] = df['Alleles_reported_Pathogenic_Likely_pathogenic'].replace(
        '-', 0)
    df['Alleles_reported_Pathogenic_Likely_pathogenic'] = df['Alleles_reported_Pathogenic_Likely_pathogenic'].astype(
        float)
    df_capped = df[df['Alleles_reported_Pathogenic_Likely_pathogenic'] < capped]
    # Plot histogram of 'Alleles_reported_Pathogenic_Likely_pathogenic' column
    plt.hist(df_capped['Alleles_reported_Pathogenic_Likely_pathogenic'].dropna(), bins=30, alpha=0.5)
    plt.title('Alleles Reported Pathogenic / Likely Pathogenic in ClinVar by gene, capped at ' + str(capped))
    plt.xlabel('Number of alleles')
    plt.ylabel('Number of genes')
    plt.grid(True)
    if SAVING_PLOT:
        plt.savefig(FIG_PATH + 'plt_hist_clinvar_patho.png')
    else:
        plt.show()


def plt_scatter_clinvar_esshap1(df_genes):
    # Create scatter plot
    plt.figure(figsize=(10, 6))
    plt.scatter(df_genes['ratio'], df_genes['Alleles_reported_Pathogenic_Likely_pathogenic'], color='gray',
                label='Genes filtered out', alpha=0.3)

    df_genes_int = df_genes[
        (df_genes['ratio'] <= 0.25) & (df_genes['Alleles_reported_Pathogenic_Likely_pathogenic'] >= 10)]
    plt.scatter(df_genes_int['ratio'], df_genes_int['Alleles_reported_Pathogenic_Likely_pathogenic'], color='blue',
                label='Genes intermediate')

    df_genes_stringent = df_genes[
        (df_genes['ratio'] <= 0.2) & (df_genes['Alleles_reported_Pathogenic_Likely_pathogenic'] >= 20)]
    plt.scatter(df_genes_stringent['ratio'], df_genes_stringent['Alleles_reported_Pathogenic_Likely_pathogenic'],
                color='green', label='Genes kept')

    # Annotate gene names for specific points
    for index, row in df_genes_int.iterrows():
        plt.annotate(row['GENE_SYMBOL'], (row['ratio'], row['Alleles_reported_Pathogenic_Likely_pathogenic']),
                     textcoords="offset points", xytext=(0, 10), ha='center', fontsize=8, color='black')

    # Add dashed lines for different filtering options
    plt.axvline(x=0.2, color='green', linestyle='--', label='Hap1 essentiality ratio < 0.2')
    plt.axvline(x=0.25, color='blue', linestyle='--', label='Hap1 essentiality ratio < 0.25')
    plt.axhline(y=10, color='lightblue', linestyle='--', label='Likely / pathogenic variant in ClinVar > 10')
    plt.axhline(y=15, color='blue', linestyle='--', label='Likely / pathogenic variant in ClinVar > 15')
    plt.axhline(y=20, color='green', linestyle='--', label='Likely / pathogenic variant in ClinVar > 20')

    # Display the number of genes left after each filter combination
    n_0210 = \
        df_genes[(df_genes['ratio'] <= 0.2) & (df_genes['Alleles_reported_Pathogenic_Likely_pathogenic'] >= 10)].shape[
            0]
    n_0215 = \
        df_genes[(df_genes['ratio'] <= 0.2) & (df_genes['Alleles_reported_Pathogenic_Likely_pathogenic'] >= 15)].shape[
            0]
    n_0220 = \
        df_genes[(df_genes['ratio'] <= 0.2) & (df_genes['Alleles_reported_Pathogenic_Likely_pathogenic'] >= 20)].shape[
            0]

    n_02510 = \
        df_genes[(df_genes['ratio'] <= 0.25) & (df_genes['Alleles_reported_Pathogenic_Likely_pathogenic'] >= 10)].shape[
            0]
    n_02515 = \
        df_genes[(df_genes['ratio'] <= 0.25) & (df_genes['Alleles_reported_Pathogenic_Likely_pathogenic'] >= 15)].shape[
            0]
    n_02520 = \
        df_genes[(df_genes['ratio'] <= 0.25) & (df_genes['Alleles_reported_Pathogenic_Likely_pathogenic'] >= 20)].shape[
            0]

    plt.text(0.2, 10, f'n= {n_0210}', fontsize=12, color='black', ha='right', va='bottom')
    plt.text(0.25, 10, f'n={n_02510}', fontsize=12, color='black', ha='right', va='bottom')
    plt.text(0.2, 15, f'n={n_0215}', fontsize=12, color='black', ha='right', va='bottom')
    plt.text(0.25, 15, f'n={n_02515}', fontsize=12, color='black', ha='right', va='bottom')
    plt.text(0.2, 20, f'n={n_0220}', fontsize=12, color='black', ha='right', va='bottom')
    plt.text(0.25, 20, f'n={n_02520}', fontsize=12, color='black', ha='right', va='bottom')

    plt.title('Number of ClinVar Pathogenic vs Essentiality Ratio: filtering options (N genes =' + str(
        df_genes.shape[0]) + ')')
    plt.xlabel('Essentiality Ratio')
    plt.ylabel('Number of ClinVar Pathogenic')
    plt.legend()
    plt.grid(True)
    if SAVING_PLOT:
        plt.savefig(FIG_PATH + 'plt_scatter_clinvar_esshap1.png')
    else:
        plt.show()


def plt_pxscatter_clinvar_esshap1(df_genes):
    fig = px.scatter(df_genes, x='ratio', y='Alleles_reported_Pathogenic_Likely_pathogenic', color='CEGs',
                     custom_data=['GENE_SYMBOL', 'ENSEMBL_ID', 'ratio', 'higher_gel_status', 'lof.oe', 'mis.oe',
                                  'syn.oe', 'CEGs', 'CSEGs'])

    # Add dashed lines for different filtering options
    fig.add_hline(y=10, line_dash="dash", line_color="lightblue",
                  annotation_text="Likely / pathogenic variant in ClinVar > 10")
    fig.add_hline(y=15, line_dash="dash", line_color="blue",
                  annotation_text="Likely / pathogenic variant in ClinVar > 15")
    fig.add_hline(y=20, line_dash="dash", line_color="green",
                  annotation_text="Likely / pathogenic variant in ClinVar > 20")
    fig.add_vline(x=0.2, line_dash="dash", line_color="green", annotation_text="Hap1 essentiality ratio < 0.2")
    fig.add_vline(x=0.25, line_dash="dash", line_color="blue", annotation_text="Hap1 essentiality ratio < 0.25")

    # Display the number of genes left after each filter combination
    n_0210 = \
        df_genes[(df_genes['ratio'] <= 0.2) & (df_genes['Alleles_reported_Pathogenic_Likely_pathogenic'] >= 10)].shape[
            0]
    n_0215 = \
        df_genes[(df_genes['ratio'] <= 0.2) & (df_genes['Alleles_reported_Pathogenic_Likely_pathogenic'] >= 15)].shape[
            0]
    n_0220 = \
        df_genes[(df_genes['ratio'] <= 0.2) & (df_genes['Alleles_reported_Pathogenic_Likely_pathogenic'] >= 20)].shape[
            0]

    n_02510 = \
        df_genes[(df_genes['ratio'] <= 0.25) & (df_genes['Alleles_reported_Pathogenic_Likely_pathogenic'] >= 10)].shape[
            0]
    n_02515 = \
        df_genes[(df_genes['ratio'] <= 0.25) & (df_genes['Alleles_reported_Pathogenic_Likely_pathogenic'] >= 15)].shape[
            0]
    n_02520 = \
        df_genes[(df_genes['ratio'] <= 0.25) & (df_genes['Alleles_reported_Pathogenic_Likely_pathogenic'] >= 20)].shape[
            0]

    # Annotate the number of genes left after each filter combination
    fig.add_annotation(x=0.2, y=10, text=f'n= {n_0210}', showarrow=False)
    fig.add_annotation(x=0.25, y=10, text=f'n={n_02510}', showarrow=False)
    fig.add_annotation(x=0.2, y=15, text=f'n={n_0215}', showarrow=False)
    fig.add_annotation(x=0.25, y=15, text=f'n={n_02515}', showarrow=False)
    fig.add_annotation(x=0.2, y=20, text=f'n={n_0220}', showarrow=False)
    fig.add_annotation(x=0.25, y=20, text=f'n={n_02520}', showarrow=False)

    fig.update_layout(title='Scatter Plot of Essentiality Ratio vs. Number of ClinVar Pathogenic',
                      xaxis_title='Essentiality Ratio',
                      yaxis_title='Number of ClinVar Pathogenic',
                      legend_title='Core Essential Gene',
                      hovermode='closest',
                      showlegend=True,
                      template='plotly_white')
    fig.update_traces(
        hovertemplate="<br>".join([
            "%{customdata[0]}",
            "ENSEMBL ID: %{customdata[1]}",
            "HAP1 essentiality ratio: %{customdata[2]}",
            "Higher GEL status: %{customdata[3]}",
            "O/E LOF: %{customdata[4]}",
            "O/E missense: %{customdata[5]}",
            "O/E synonymous: %{customdata[6]}",
            "Cores Essential Gene: %{customdata[7]}",
            "Cancer Specific EG: %{customdata[8]}",
        ])
    )

    if SAVING_PLOT:
        fig.write_html(FIG_PATH + 'px_filtering_clinvar_ess_ratio.html')
    else:
        fig.show()


def plt_scatter_clinvar_lofoe(df_genes, threshold_ess):
    # Create scatter plot
    plt.figure(figsize=(12, 12))
    scatter = plt.scatter(df_genes['lof.oe'], df_genes['ratio'], c=df_genes['CEGs'], alpha=1)

    # Annotate gene names for specific points
    for index, row in df_genes.iterrows():
        plt.annotate(row['GENE_SYMBOL'], (row['lof.oe'], row['ratio']),
                     textcoords="offset points", xytext=(0, 10), ha='center', fontsize=8, color='black')
    # Color bar legend
    cbar = plt.colorbar(scatter)
    cbar.set_label('Core essential genes')
    plt.title('HAP1 Essentiality ratio vs pLOF o/e: on selected gene (N genes =' + str(
        df_genes.shape[0]) + ')\nHAP1 essential ratio <= ' + str(threshold_ess))
    plt.xlabel('o/e pLOF')
    plt.ylabel('HAP1 Essentiality ratio')
    plt.legend()
    plt.grid(True)

    if SAVING_PLOT:
        plt.savefig(FIG_PATH + 'plt_scatter_clinvar_lofoe.png')
    else:
        plt.show()


def plt_pxscatter_clinvar_lofoe(df_genes):
    # Convert dictionary objects to JSON strings
    df_genes['mode_of_pathogenicity'] = df_genes['mode_of_pathogenicity'].apply(json.dumps)
    df_genes['mode_of_inheritance'] = df_genes['mode_of_inheritance'].apply(json.dumps)

    fig = px.scatter(df_genes, x='lof.oe', y='ratio', color='CEGs',
                     custom_data=['GENE_SYMBOL', 'ENSEMBL_ID', 'ratio', 'higher_gel_status', 'lof.oe', 'mis.oe',
                                  'syn.oe', 'CEGs', 'CSEGs', 'mode_of_pathogenicity', 'mode_of_inheritance',
                                  'Alleles_reported_Pathogenic_Likely_pathogenic'])
    # Add annotations for GENE_SYMBOL next to the points
    for i, row in df_genes.iterrows():
        fig.add_annotation(x=row['lof.oe'], y=row['ratio'] + 0.003,
                           text=row['GENE_SYMBOL'], showarrow=False, font=dict(size=10))

    # Display the number of genes left after each filter combination
    n_05 = df_genes[(df_genes['lof.oe'] <= 0.5)].shape[0]
    n_06 = df_genes[(df_genes['lof.oe'] <= 0.6)].shape[0]
    n_07 = df_genes[(df_genes['lof.oe'] <= 0.7)].shape[0]
    n_08 = df_genes[(df_genes['lof.oe'] <= 0.8)].shape[0]
    n_09 = df_genes[(df_genes['lof.oe'] <= 0.9)].shape[0]

    # Annotate the number of genes left after each filter combination
    # Add dashed lines for different filtering options
    fig.add_vline(x=0.5, line_dash="dash", line_color="black", annotation_text=f'n= {n_05}')
    fig.add_vline(x=0.7, line_dash="dash", line_color="black", annotation_text=f'n= {n_07}')
    fig.add_vline(x=0.6, line_dash="dash", line_color="black", annotation_text=f'n= {n_06}')
    fig.add_vline(x=0.8, line_dash="dash", line_color="black", annotation_text=f'n= {n_08}')
    fig.add_vline(x=0.9, line_dash="dash", line_color="black", annotation_text=f'n= {n_09}')

    fig.update_layout(title='Scatter Plot of HAP1 Essentiality ratio vs O/E pLOF (N=' + str(df_genes.shape[0]) + ')',
                      xaxis_title='O/E pLOF',
                      yaxis_title='HAP1 Essentiality ratio',
                      legend_title='Core Essential Gene',
                      hovermode='closest',
                      showlegend=True,
                      template='plotly_white')
    fig.update_traces(
        hovertemplate="<br>".join([
            "%{customdata[0]}",
            "ENSEMBL ID: %{customdata[1]}",
            "HAP1 essentiality ratio: %{customdata[2]}",
            "Higher GEL status: %{customdata[3]}",
            "O/E LOF: %{customdata[4]}",
            "O/E missense: %{customdata[5]}",
            "O/E synonymous: %{customdata[6]}",
            "Cores Essential Gene: %{customdata[7]}",
            "Cancer Specific EG: %{customdata[8]}",
            "mode_of_pathogenicity: %{customdata[9]}",
            "mode_of_inheritance: %{customdata[10]}",
            "ClinVar Pathogenic Likely pathogenic: %{customdata[11]}",
        ]))
    if SAVING_PLOT:
        fig.write_html(FIG_PATH + 'px_filtering_clinvar_oe_plof.html')
    else:
        fig.show()


# MAIN ------------------------------------------------------------------------------------------------------------
# Create the log file
if SAVING_DF:
    log_file_path = 'data/out_df/v3_' + DATE + '_gene_preprocess.log'
    log_file = open(log_file_path, 'a')
    sys.stdout = log_file

# preprocessing
panel_df = cleaning_panel_bed(panel_df)


panel_df.rename(columns={'ensemblGenes': 'ENSEMBL_ID', 'geneSymbol': 'geneSymbol_panel_app'}, inplace=True)

oe_df = cleaning_oe(df_oe_all)
# print_info_gene_panel(panel_df['geneSymbol'].unique(), panel_df, False)
HAP1_blomen_ess_df = HAP1_blomen_df[HAP1_blomen_df.selected == "YES"]

# Intersection between HAP1 and GenePanel
print('\nNumber of essential genes in HAP1 Blomen:', HAP1_blomen_ess_df.shape[0])
intersection_symbol = set(HAP1_blomen_ess_df['GENE_SYMBOL']).intersection(panel_df['geneSymbol_panel_app'].unique())
intersection_ensbl = set(HAP1_blomen_ess_df['ENSEMBL_ID']).intersection(panel_df['ENSEMBL_ID'].unique())
print('Number of genes HAP1 essential and in gene panel, match by gene symbol', len(intersection_symbol))
print('Number of genes HAP1 essential and in gene panel, match by gene ensemble', len(intersection_ensbl))

# Create df of the 865 genes of interest ---------------------------------------------------------------------------
df_genes = HAP1_blomen_ess_df[HAP1_blomen_ess_df['ENSEMBL_ID'].isin(intersection_ensbl) | (
    HAP1_blomen_ess_df['GENE_SYMBOL'].isin(intersection_symbol))]
print('Tot number of genes HAP1 essential and in gene panel, matching by gene symbol or by ENSMB', df_genes.shape[0])

# Updating the DataFrame
df_genes.loc[df_genes['ENSEMBL_ID'] == 'ENSG00000131795', 'GENE_SYMBOL'] = 'RBM8A'
df_genes.loc[df_genes['ENSEMBL_ID'] == 'ENSG00000131795', 'ENSEMBL_ID'] = 'ENSG00000265241'

df_genes.loc[df_genes['ENSEMBL_ID'] == 'ENSG00000108278', 'GENE_SYMBOL'] = 'ZNHIT3'
df_genes.loc[df_genes['ENSEMBL_ID'] == 'ENSG00000108278', 'ENSEMBL_ID'] = 'ENSG00000273611'
updated_inters_ensbl = list(intersection_ensbl) + ['ENSG00000265241', 'ENSG00000273611']
print(updated_inters_ensbl)
df_genes['higher_gel_status'] = df_genes['ENSEMBL_ID'].map(get_gene_higher_gel_status(updated_inters_ensbl, panel_df))

panel_df_unique = panel_df.drop_duplicates('ENSEMBL_ID', keep='last')
df_genes = df_genes.merge(panel_df_unique[['ENSEMBL_ID', 'geneSymbol_panel_app']], on='ENSEMBL_ID', how='left')
print('-- df_genes with intersection Hap1 panelApp done')
#print_info_gene_panel(df_genes['ENSEMBL_ID'].unique(), panel_df, False)

# Add POSITION start end in hg38
gene_gtf_df.rename(columns={'gene_id': 'ENSEMBL_ID'}, inplace=True)
df_genes = df_genes.merge(gene_gtf_df, on='ENSEMBL_ID', how='left')
nan_count_limit = df_genes['start'].isna().sum()
print('-- df_genes adding hg38 position done')
print("Number of NaN values in the 'start' column:", nan_count_limit, 'over', df_genes.shape[0])

# add GNOMAD 'lof.oe', 'mis.oe', 'syn.oe'
df_genes['GENE_SYMBOL'] = df_genes['gene_name']
df_genes = df_genes.merge(oe_df[['GENE_SYMBOL', 'lof.oe', 'mis.oe', 'syn.oe']], on='GENE_SYMBOL', how='left')
print('-- df_genes adding Gnomad info done')
nan_count_limit = df_genes['mis.oe'].isna().sum()
print("Number of NaN values in the 'mis.oe' column:", nan_count_limit, 'over', df_genes.shape[0])

# Adding OGEE essentiality
ogee_CSEGs_ens_list, ogee_CEGs_ens_list, ogee_CSEGs_gene_list, ogee_CEGs_gene_list = get_ensembl_list(df_ogee)
df_genes['CSEGs'] = (
        df_genes['ENSEMBL_ID'].isin(ogee_CSEGs_ens_list) | df_genes['GENE_SYMBOL'].isin(ogee_CSEGs_gene_list))
df_genes['CEGs'] = (df_genes['ENSEMBL_ID'].isin(ogee_CEGs_ens_list) | df_genes['GENE_SYMBOL'].isin(ogee_CEGs_gene_list))
print('-- df_genes adding CSEGs and CSEGs info done')

# Adding ClinVar summary
clinvar_gene_specific_summary.rename(columns={'#Symbol': 'GENE_SYMBOL'}, inplace=True)
df_genes = df_genes.merge(clinvar_gene_specific_summary, on='GENE_SYMBOL', how='left')
df_genes['Alleles_reported_Pathogenic_Likely_pathogenic'] = df_genes[
    'Alleles_reported_Pathogenic_Likely_pathogenic'].replace('-', 0)
df_genes['Alleles_reported_Pathogenic_Likely_pathogenic'] = df_genes[
    'Alleles_reported_Pathogenic_Likely_pathogenic'].astype(float)

# Manual fill NaN for ClinVar
dict_manual_matching = {'ENSG00000145375': [166378, 996, 793, 372, 104, 613940, 357, 20],  # AFG2A
                        'ENSG00000171763': [79029, 123, 93, 93, 17, 619578, 52, 2]}  # AFG2B
for ensembl_id, replacement_values in dict_manual_matching.items():
    mask = df_genes['ENSEMBL_ID'] == ensembl_id
    if mask.any():
        df_genes.loc[mask, ['GeneID', 'Total_submissions', 'Total_alleles',
                            'Submissions_reporting_this_gene', 'Alleles_reported_Pathogenic_Likely_pathogenic',
                            'Gene_MIM_number', 'Number_uncertain', 'Number_with_conflicts']] = replacement_values
print('-- df_genes adding ClinVar info done')
nan_count_limit = df_genes['Alleles_reported_Pathogenic_Likely_pathogenic'].isna().sum()
print("Number of NaN values in the 'Alleles_reported_Pathogenic_Likely_pathogenic' column:", nan_count_limit, 'over',
      df_genes.shape[0])

# Add mode of inheritance /pathogenicity from gene panels
mode_of_pathogenecity_dict, mode_of_inheritance_dict = dict_modeofinheritance_modeofpatho(panel_df)
df_genes['mode_of_pathogenicity'] = df_genes['ENSEMBL_ID'].map(mode_of_pathogenecity_dict)
df_genes['mode_of_inheritance'] = df_genes['ENSEMBL_ID'].map(mode_of_inheritance_dict)
print('-- df_genes adding mode_of_pathogenicity and mode_of_inheritance info done')

# PLOT + Selection ----------------------------------------------------------------------------------------------------
if PLOT:
    plt_swarm_essratio_highergelstatus(df_genes)
    plt_scatter_oe_bhap1ratio(df_genes)
    plt_hist_oe_by_ogeev3_ess(df_genes)
    plt_hist_ratio_by_ogeev3_ess(df_genes)
    plt_hist_clinvar_path(df_genes)

# Remove
print('\n\nSELECTION -----')

print('-- Initial number of genes:', len(df_genes['ENSEMBL_ID']))
df_genes_no_selection = df_genes.copy()
# keep only level confidence ==3
df_genes = df_genes[(df_genes.higher_gel_status == 3)]
print('-- genes with higher GEL status = 1,2 are removed panel red, amber dropped')
print('-- Number of genes in current list:', len(df_genes['ENSEMBL_ID']))

if PLOT:
    plt_pxscatter_clinvar_esshap1(df_genes)
    plt_scatter_clinvar_esshap1(df_genes)

# Print NaN values
print('Total number of NaN for gene selected =', df_genes.isna().sum().sum())
df_panel_subset = print_info_gene_panel(df_genes['ENSEMBL_ID'].unique(), panel_df, False)
# merge df_panel and df_genes

print('\nfinal df_genes_no_selection : ')
print(df_genes.head())
print(df_genes.shape)
print(df_panel_subset.head())
print(df_panel_subset.shape)
print(' Number of unique genes:', len(df_panel_subset['geneSymbol_panel_app'].unique()))
print(' Number of unique genes:', len(df_panel_subset['ENSEMBL_ID'].unique()))
# Merge the two DataFrames based on ENSEMBL_ID
merged_df = pd.merge(df_panel_subset, df_genes, on='ENSEMBL_ID', how='inner')
print(merged_df.head())
print(merged_df.shape)
path_merged_df_df = '/Users/terwagc/PycharmProjects/EG_visualisation/data/output_df/' + DATE + '_' + str(
    merged_df.shape[0]) + '_genes_disease_selected.tsv'
merged_df.to_csv(path_merged_df_df, sep='\t')

if SAVING_DF:
    path_subset_df = 'data/out_df/' + DATE + '_' + str(
        df_genes.shape[0]) + '_genes_selected.tsv'
    path_all_df = 'data/out_df/' + DATE + '_' + str(
        df_genes_no_selection.shape[0]) + '_genes_selected.tsv'
    df_genes.to_csv(path_subset_df, sep='\t')
    df_genes_no_selection.to_csv(path_all_df, sep='\t')
    print("\nSaving subset df in ", path_subset_df)
    print("Saving Full df in ", path_all_df)
if SAVING_DF:
    log_file.close()
