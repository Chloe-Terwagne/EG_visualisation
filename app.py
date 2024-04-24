"""
    Project :  variants intersection visualisaiton platform for BRCA1
    Description: Show interactively variant present, in both, raw BRCA1 SGE data from G. Findlay et al., 2018 and the UKB allele table
    Name : Chloé Terwagne
    date : 25 April 2023
    Python version 3.10
"""
# IMPORT ---------------------------------------------------------------

import numpy as np
from dash import Dash, dcc, html, Output, Input
import dash_bootstrap_components as dbc
import plotly.express as px
import pandas as pd
from protein_folding import create_style_3d
import dash_bio as dashbio
from dash_bio.utils import PdbParser
from dash.development.base_component import Component, _explicitize_args

# FONT & COLOR  ---------------------------------------------------------------
font_list = ["Arial", "Balto", "Courier New", "Droid Sans", "Droid Serif", "Droid Sans Mono", "Gravitas One",
             "Old Standard TT", "Open Sans", "Overpass", "PT Sans Narrow", "Raleway", "Times New Roman"]
idx_font = 0  # UCL font
# Color Model
pink_vibrant = 'rgb(172, 20, 90)'  # UCL color
blue_vibrant = 'rgb(52, 198, 198)'  # UCL color
yellow_vibrant = 'rgb(255, 202, 56)'  # UCL color

yellow_mutated = 'rgb(218, 214, 202)'  # UCL color
pink_mutated = 'rgb(222, 184, 195)'  # UCL color
purple_mutated = "rgb(198,176,188)"  # UCL color


# pd display
templates = 'plotly_dark'
pd.set_option('display.width', 900)
pd.set_option('display.max_columns', 200)
pd.set_option("display.max_rows", None)


# FUNCTION  -----------------------------------------------------------------------------------------------------------

# MAIN ---------------------------------------------------------------------------------------------------------------

# cleaning data

df = pd.read_csv(
    "https://github.com/Chloe-Terwagne/VarEffectViz/blob/main/df/merged_brca1_sge_ukb_2023_04_21.csv?raw=true")

# glossary padding
cell_style = {'padding-bottom': '20px', 'font-weight': 'bold', 'color': yel}  # more title type
bottom_style = {'padding-bottom': '15px'}  # more cell type

# Create abbreviation hover
sge_abbr = html.Abbr('SGE', title=' Saturation Genome Editing ')
ukb_abbr = html.Abbr('UKB', title=' UKBiobank ')
cadd_abbr = html.Abbr('CADD', title=' Combined Annotation Dependent Depletion ')
ac_abbr = html.Abbr('AC', title=' Allele Count ')

# Launch app------------------------------------------------------------------------------------------------
app = Dash(__name__, external_stylesheets=[dbc.themes.DARKLY], suppress_callback_exceptions=True,
           meta_tags=[{'name': 'viewport', 'content': 'width=device-width, initial-scale=1.0'}])
server = app.server

# Build your components------------------------------------------------------------------------------------------------
# 3D parsing & styling
parser = PdbParser('https://github.com/Chloe-Terwagne/VarEffectViz/blob/main/df/AF-P38398-F1-model_v4.pdb?raw=true')
# from https://alphafold.ebi.ac.uk/entry/P38398
data = parser.mol3d_data()
styles = create_style_3d(
    df, 'minmax_neg_func_score', data['atoms'], visualization_type='cartoon', color_element='residue_score')
brca1_3D = dashbio.Molecule3dViewer(id='dashbio-default-molecule3d',modelData=data,styles=styles,backgroundOpacity=0,
    selectionType='residue',backgroundColor="black",height=500, width=1250, zoom=dict(factor=1.9,
                                                                    animationDuration=30000, fixedPath=False))
overview_title = dcc.Markdown(children='', style=dict(font_family=font_list[idx_font], font_color=yel))
overview_display = dcc.RadioItems(options=["Variants aggregated by position", "Variants expanded by nucleotide type"],
                                  value='Variants aggregated by position', labelClassName="custom-text p-3")
items = ['Clinvar', 'Consequence', "SGE", "UKB"]
overview_dropdown = dcc.Dropdown(options=['Clinvar', 'Consequence', "SGE", "UKB"],
                                 value='Consequence', clearable=False, className='my-custom-dropdown')
overview_graph = dcc.Graph(figure={}, config={'staticPlot': False, 'scrollZoom': False, 'doubleClick': 'reset',
                                              'showTips': True,'displayModeBar': 'hover', 'displaylogo': False,
                                              'modeBarButtonsToRemove': ['lasso2d', 'zoomIn2d', 'zoomOut2d',
                                                                         'autoScale2d'], }, selectedData=None)
color_blind_option = BooleanSwitch(on=False, size=30,
                                   label=dict(label="color blind friendly", style=dict(font_color=yel)),
                                   color=mid_purple, labelPosition="left")
exon_option = BooleanSwitch(on=False, size=30, label=dict(label="exon highlighting", style=dict(font_color=yel)),
                            color=yellow, labelPosition="left")
three_d_graph = dcc.Graph(figure={},
                          config={'staticPlot': False, 'scrollZoom': True, 'doubleClick': 'reset', 'showTips': True,
                                  'displayModeBar': False, 'watermark': False})
three_d_title = dcc.Markdown(children='all variant')
clinvar_hist_graph_sge = dcc.Graph(figure={}, config={'staticPlot': False, 'scrollZoom': False, 'doubleClick': 'reset',
                                                      'showTips': True, 'displayModeBar': False, 'watermark': False})
clinvar_hist_graph_ac = dcc.Graph(figure={}, config={'staticPlot': False, 'scrollZoom': False, 'doubleClick': 'reset',
                                                     'showTips': True, 'displayModeBar': False, 'watermark': False})
clinvar_hist_graph_cadd = dcc.Graph(figure={}, config={'staticPlot': False, 'scrollZoom': False, 'doubleClick': 'reset',
                                                       'showTips': True, 'displayModeBar': False, 'watermark': False,
                                                       })
mol_viewer_colorbar = dcc.Graph(figure={}, config={'staticPlot': True, 'scrollZoom': False, 'showTips': False,
                                                   'displayModeBar': False, 'watermark': False})

github_link = html.Div([
    html.A(
        id='gh-link',
        children=['View on GitHub'],
        href="https://github.com/Chloe-Terwagne/VarEffectViz",
        style={'color': yel, 'border': yel, "text-decoration": 'none'},
        target="_blank"
    ),
    html.Img(src='https://github.com/Chloe-Terwagne/VarEffectViz/blob/main/src/assets/GitHub-Mark-64px.png?raw=true',
             style={'height': '50px', 'margin-left': '10px'},
             ),
], style={
    'background': 'black',
    'color': 'white',
    'height': '80px',
    'width': '100%',
    'border': 'solid 2px white',
    'display': 'flex',
    'align-items': 'center',
    'justify-content': 'center',
    'padding': '20px',
    "margin-top": "180px",
    "margin-bottom": "10px"
})

text_abreviation = dbc.Card(
    [
        dbc.CardImg(
            src="https://github.com/Chloe-Terwagne/VarEffectViz/blob/main/src/assets/ucl-banner-port-stone-rgb-lg.png?raw=true",
            top=True),
        dbc.CardBody(
            [
                html.H4("Quick Resources", className="app-controls-block",
                        style={"font-family": "Garamond", 'margin-top': '10px', 'margin-bottom': '20px',
                               'font-size': '18pt'}),
                html.P(
                    "Hover over the acronym to reveal the full name.",
                    className="card-text", style={'font-size': '9pt'}
                ),

                html.Div([
                    sge_abbr, html.Br(),
                    ukb_abbr, html.Br()], style={'float': 'left', 'width': '50%'}),

                html.Div([
                    cadd_abbr, html.Br(),
                    ac_abbr, html.Br()], style={'float': 'left', 'width': '50%'}),

                github_link,
                dbc.CardLink(["G.Findlay", html.Em(" et al."), ', 2018'],
                             href="https://www.nature.com/articles/s41586-018-0461-z",
                             target="_blank", className='custom-link'),
                dbc.CardLink("Genome Function Lab", href="https://www.crick.ac.uk/research/labs/greg-findlay/",
                             target="_blank", className='custom-link'),
                dbc.CardLink("UKB initiative",
                             href="https://www.ukbiobank.ac.uk/learn-more-about-uk-biobank/about-us/", target="_blank",
                             className='custom-link'),
                # dbc.CardLink("AlphaFold",
                #              href="https://alphafold.ebi.ac.uk/entry/P38398", target="_blank",className='custom-link'),#https://www.ncbi.nlm.nih.gov/clinvar/intro/
                # dbc.CardLink("CADD",
                #              href="https://cadd.gs.washington.edu/info", target="_blank",className='custom-link'),
                # dbc.CardLink("Clinvar",
                #              href="https://www.ncbi.nlm.nih.gov/clinvar/intro/", target="_blank",className='custom-link'),
            ], className='my-custom-background'
        )
    ],
    style={"height": "500px", 'background-color': 'rgba(41,41,41,0)', 'border': 'solid 2px rgb(214,210,196)',
           'border-radius': '20px',
           'overflow': 'hidden'}
)

# Customize Layout--------------------------------------------------------------------------------------------
row_style = {'display': 'flex', 'flex-wrap': 'wrap', 'align-items': 'stretch'}
app.layout = \
    dbc.Container([
        dbc.Row([html.Br()]),

        dbc.Row([
            dbc.Col(html.H1("Variant effect in BRCA1 gene", className='custom-h1'), width={'size': 7, 'offset': 2}, ),
            dbc.Col([
                dbc.Row(color_blind_option, className="my-custom-switch"),
                dbc.Row(exon_option, className="my-custom-switch")], width={'size': 2}, align='right')
        ], justify='between'),
        dbc.Row([html.Br()]),
        dbc.Row([html.Br()]),
        dbc.Row([dbc.Col(overview_display, width={'size': 6}),
                 dbc.Col([overview_dropdown], width={'size': 2}, align='right', className='my-custom-dropdown'),
                 ], justify='between'),
        dbc.Row([
            dbc.Col(overview_graph, width=12)
        ], justify='around'),
        dbc.Row([html.Br()]),
        dbc.Row([
            dbc.Col([
                dbc.Card(
                    [
                        html.Div(
                            id='body',
                            className='app-body',
                            children=[
                                html.Div(
                                    id='desc-control-tabs',
                                    className='control-tabs',
                                    children=[
                                        dcc.Tabs(id='about-tabs', value='what-is', children=[
                                            dcc.Tab(
                                                label='About',
                                                value='what-is',
                                                children=html.Div(className='control-tab', children=[
                                                    html.H4(className='app-controls-block',
                                                            children='What is VarEffectViz?'),
                                                    html.P(
                                                        'VarEffectViz is a visualizer that presents variant effect interpretation for the BRCA1 gene using multiple sources of data, including allele count in UKBiobank, the CADD computational score, AlphaFold  protein structure prediction and the results of saturation genome editing experiments. '),
                                                    html.P(
                                                        'By integrating data from these diverse sources, VarEffectViz provides a more comprehensive view of the pathogenicity of each variant. This approach allows researchers and clinicians to better understand the impact of BRCA1 genetic variants, which has important implications for cancer risk and prevention.'),
                                                    html.P([
                                                               'The "Glossary" tab explains key concepts related to genetic variants and their effects. ',
                                                               html.A("A two-minute demo video",
                                                                      href="https://youtu.be/t9ady9CxtI0",
                                                                      target="_blank",
                                                                      style={"font-family": "Garamond", "color": yel,
                                                                             "text-decoration": "none"}
                                                                      ),
                                                               ' showcasing the multiple functionalities of this visualization board is also available. Please note that since the board is deployed as a free version, some updates may take up to 10-15 seconds.'])
                                                ])
                                            ),
                                            dcc.Tab(
                                                label='Variant effect',
                                                value='var-effect',
                                                children=html.Div(className='control-tab', children=[
                                                    html.H4(className='app-controls-block',
                                                            children='The variant interpretation challenges'),
                                                    html.P(
                                                        'Every human being has a unique genetic code that is responsible for many of their physical and biological characteristics. The genetic code is made up of DNA, which is organized into distinct units called genes. When changes occur in the DNA sequence, these changes are known as variants. Some variants have no effect on health or function, while others can lead to disease or altered biological processes.'),
                                                    html.P(
                                                        'To better understand the potential effects of a variant, scientists and clinicians use a process called variant interpretation. This involves analyzing the genetic changes to determine whether they are benign or pathogenic. A combination of different approaches is used in this process, including assessing the frequency of the variant in the general population, using computational algorithms to predict the effect of the variant on biological function, and performing experimental studies to validate these predictions.'),
                                                    html.P(
                                                        'However, variant interpretation is not a perfect process, and each approach has its own strengths and limitations. Therefore, to achieve the best possible understanding of genetic variation, researchers and clinicians must combine different approaches to obtain a more comprehensive view of the potential effects of a variant.'),
                                                    html.P(
                                                        'By using multiple sources of data, variant interpretation can provide valuable insights into the role of genetics in health and disease. It can help guide medical diagnosis, treatment, and prevention by identifying genetic changes that may increase the risk of disease or influence response to treatment. Overall, variant interpretation is a crucial tool in the field of genetics and is essential for advancing our understanding of the role of genetics in health and disease.'),
                                                ])
                                            ),
                                            dcc.Tab(
                                                label='Data',
                                                value='data-resource',
                                                children=html.Div(className='control-tab', children=[
                                                    html.H4(className='app-controls-block', children='Data sources'),
                                                    html.P([
                                                        "The variants featured in this visualisation board have been observed at least once in the UKBiobank cohort and have undergone experimental testing through saturation genome editing, as described in the ",
                                                        html.Em("Nature"),
                                                        " publication ",
                                                        html.A(
                                                            "Accurate classification of BRCA1 variants with saturation genome editing",
                                                            href="https://www.nature.com/articles/s41586-018-0461-z",
                                                            target="_blank",
                                                            style={"font-family": "Garamond", "color": yel,
                                                                   "text-decoration": "none"}
                                                        ),
                                                        " by Greg Findlay ",
                                                        html.Em("et al."),
                                                        " in 2018. Further information regarding the data's origin and implications can be found below."
                                                    ]),
                                                    html.P(
                                                        "Saturation Genome Editing is a lab technique for testing the effects of genetic variants on protein function in cells. It involves systematically introducing mutations into the DNA sequence of a gene using a genome editing tool such as CRISPR-Cas9, and then assessing the resulting changes in the function of the protein that the gene encodes."),
                                                    html.P(
                                                        "Saturation Genome Editing function score is assigned to a genetic variant based on its effect on cell survival rates when the variant is introduced into a cell using the Saturation Genome Editing technique. The score reflects the degree to which the variant disrupts the normal function of the protein encoded by the gene. Saturation Genome Editing function scores are experimental measures of variant pathogenicity."),
                                                    html.P(
                                                        "UK Biobank is the biggest human sequencing project to date, containing genetic and health-related data from over 500,000 participants in the United Kingdom. UK Biobank is a valuable resource for researchers studying the genetic basis of diseases, as it provides a large sample size and diverse set of genetic and health-related data."),

                                                    html.P(
                                                        "A deep learning system developed by Google's DeepMind, AlphaFold, that predicts the 3D structure of a protein based on its amino acid sequence. The predicted structure can provide valuable information for variant interpretation, as variants can disrupt protein folding and binding, potentially affecting the function of the protein."),
                                                    html.P(
                                                        "CADD is a tool that predicts the deleteriousness of genetic variants based on their similarity to known pathogenic and benign variants. The score takes into account various genomic annotations, such as conservation, functional genomics, and regulatory information, to provide a single score that reflects the likelihood of a variant being deleterious. The higher the CADD score, the more pathogenic the variant is predicted to be."),
                                                    html.P(
                                                        "ClinVar is a public database of genetic variants and their clinical significance. ClinVar is considered a highly trusted, \"gold standard\" resource for variant interpretation, as it collects and curates variant data from multiple sources, including clinical laboratories, research studies, and expert panels. However, it's important to note that the majority of genetic variants have not yet been annotated in ClinVar, meaning that they are classified as \"unknown significance\" or \"absent from ClinVar\"."),

                                                ])
                                            ),

                                            dcc.Tab(
                                                label='Glossary',
                                                value='interpe',
                                                children=html.Div(className='control-tab', children=[
                                                    html.Br(),
                                                    html.Tr([html.Td('Term', style=cell_style), html.Td(
                                                        'Definition', style=cell_style)]),
                                                    html.Tr([html.Td('BRCA1'), html.Td(
                                                        'A gene that encodes a protein involved in DNA repair and maintenance of genomic stability. Mutations in the BRCA1 gene are associated with an increased risk of developing breast and ovarian cancer.',
                                                        style=bottom_style)]),
                                                    html.Tr([html.Td('Variant'), html.Td(
                                                        'A genetic variation that occurs in an individual\'s DNA sequence. Variants can be benign or pathogenic and can affect various biological processes in the body.',
                                                        style=bottom_style)]),
                                                    html.Tr([html.Td('DNA'), html.Td(
                                                        'The molecule that carries genetic information and is present in almost all living organisms.',
                                                        style=bottom_style)]),
                                                    html.Tr([html.Td('Nucleotide'), html.Td(
                                                        'The basic building block of DNA and RNA, consisting of a sugar molecule, a phosphate group, and one of four nitrogenous bases: adenine (A), cytosine (C), guanine (G), or thymine (T) in DNA (or uracil (U) in RNA). The sequence of nucleotides in DNA encodes genetic information, while the sequence of nucleotides in RNA is used to direct protein synthesis.',
                                                        style=bottom_style)]),
                                                    html.Tr([html.Td('RNA'), html.Td(
                                                        'A molecule that plays a critical role in the transfer of genetic information from DNA to proteins. RNA is involved in a variety of biological processes, including protein synthesis and regulation of gene expression.',
                                                        style=bottom_style)]),
                                                    html.Tr([html.Td('Amino acid'), html.Td(
                                                        'The building blocks of proteins. Amino acids are linked together to form long chains that fold into specific three-dimensional shapes to carry out various biological functions.',
                                                        style=bottom_style)]),
                                                    html.Tr([html.Td('Reference genome'), html.Td(
                                                        'A standard DNA sequence used as a reference for comparing genetic variation in different individuals. It provides a basis for identifying genetic differences that may contribute to disease or other traits.',
                                                        style=bottom_style)]),
                                                    html.Tr([html.Td('Exon'), html.Td(
                                                        'A coding region of DNA that contains the instructions for making a protein.',
                                                        style=bottom_style)]),
                                                    html.Tr([html.Td('Intron'), html.Td(
                                                        'A non-coding section of DNA that separates the coding regions (exons) of a gene. Introns are removed during the process of making RNA from DNA and do not encode proteins, but can play important regulatory roles in gene expression.',
                                                        style=bottom_style)]),
                                                    html.Tr([html.Td('Allele count'), html.Td(
                                                        'The number of copies of a particular variant in a given population.',
                                                        style=bottom_style)]),
                                                    html.Tr([html.Td('Loss of Function'), html.Td(
                                                        'A genetic variation that results in a partial or complete loss of the normal function of a protein. Such variants are often associated with pathogenicity and can lead to genetic diseases or an increased risk of disease.',
                                                        style=bottom_style)]),
                                                    html.Tr([html.Td('5\' UTR'), html.Td(
                                                        'The non-coding region at the beginning of an mRNA molecule, upstream of the start codon, that plays a role in the regulation of gene expression.',
                                                        style=bottom_style)]),
                                                    html.Tr([html.Td('Missense'), html.Td(
                                                        'A genetic variation that results in a change in the amino acid sequence of a protein.',
                                                        style=bottom_style)]),
                                                    html.Tr([html.Td('Splice acceptor'), html.Td(
                                                        'A region of DNA that signals the end of an exon and is necessary for proper splicing of the pre-mRNA.',
                                                        style=bottom_style)]),
                                                    html.Tr([html.Td("Splice donor"), html.Td(
                                                        "A region of DNA that signals the start of an exon and is necessary for proper splicing of the pre-mRNA.",
                                                        style=bottom_style)]),
                                                    html.Tr([html.Td("Splice region"), html.Td(
                                                        "A region of DNA that is important for the splicing of pre-mRNA.",
                                                        style=bottom_style)]),
                                                    html.Tr([html.Td("Start lost"), html.Td(
                                                        "A genetic variation that causes the loss of the initiation codon in the coding sequence of a gene.",
                                                        style=bottom_style)]),
                                                    html.Tr([html.Td("Stop gained"), html.Td(
                                                        "A genetic variation that creates a premature stop codon in the coding sequence of a gene.",
                                                        style=bottom_style)]),
                                                    html.Tr([html.Td("Synonymous"),
                                                             html.Td(
                                                                 "A genetic variation that does not change the amino acid sequence of a protein.",
                                                                 style=bottom_style)]),
                                                ])
                                            ),
                                        ])
                                    ],
                                    style={'overflow-y': 'auto', 'max-height': '830px'}
                                )])])], style={'background-color': 'rgba(41,41,41,0)'}, xs=12, sm=12, md=6, lg=3, xl=3),

            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        dbc.Row([
                            dbc.Col(clinvar_hist_graph_sge),
                        ]),
                        dbc.Row([
                            dbc.Col(clinvar_hist_graph_ac),
                        ]),
                        dbc.Row([
                            dbc.Col(clinvar_hist_graph_cadd),
                        ])
                    ], className='my-custom-background')
                ])
            ], className='my-custom-background', xs=12, sm=12, md=6, lg=3, xl=3),
            dbc.Col([three_d_graph], xs=12, sm=12, md=12, lg=5, xl=5)
        ], style=row_style, justify='around'),
        dbc.Row([html.Br()]),

        # row 3 ----------------------
        dbc.Row([
            # 3D protein
            dbc.Col([html.Div([brca1_3D])], xs=12, sm=12, md=6, lg=7, xl=6),
            # color bar scatter plot
            dbc.Col([html.Div([html.Br(), html.Br(), mol_viewer_colorbar, html.Hr(),
                               html.Div(id='default-molecule3d-output', style={'background-color': dark_gray_transp,
                                                                               'position': 'relative',
                                                                               'z-index': '1'
                                                                               })])], className='custom-text_left',
                    xs=12, sm=12, md=3, lg=3, xl=4),
            dbc.Col([text_abreviation], xs=12, sm=4, md=3, lg=2, xl=2),
        ], style=row_style, justify='around'),
    ], fluid=True)


def histogram(x_axis, color_blind):
    if color_blind:
        colors = {'Benign': mid_green, 'Pathogenic': mid_purple}
    else:
        colors = {'Benign': mid_green, 'Pathogenic': mid_red}

    if x_axis == 'minmax_neg_func_score':
        axlis_label = 'SGE fct score'

    elif x_axis == 'cadd_score':
        axlis_label = 'CADD score'
    else:
        axlis_label = x_axis
    df_t = df
    df_t = df_t.rename(columns={'clinvar_simple': "Clinvar high confidence", "minmax_neg_func_score": "SGE fct score",
                                'cadd_score': 'CADD score'})
    fig = px.histogram(df_t, x=axlis_label, color="Clinvar high confidence", color_discrete_map=colors, marginal="rug",
                       hover_data=["var_name", '1/AC', 'cohort_allele_count', 'CADD score', 'SGE fct score'],
                       hover_name="var_name")

    fig.update_layout(barmode='overlay')
    fig.update_traces(opacity=0.75)
    fig.update_layout(
        height=260,
        margin=dict(b=0, t=40, l=1, r=2),
        plot_bgcolor=dark_gray,
        paper_bgcolor=dark_gray,
        font_family=font_list[idx_font],
        xaxis=dict(showgrid=False, visible=True, zeroline=True, title=axlis_label),
        yaxis=dict(showgrid=False, visible=True, zeroline=True, title="variant number"),
        font_color=yel)
    if x_axis == '1/AC':
        fig.update_layout(showlegend=True,
                          legend=dict(title='Clinvar:', orientation='h', yanchor='top', y=-0.3, xanchor='left', x=0), )
    else:
        fig.update_layout(showlegend=False)
    fig.update_yaxes(showgrid=False, row=2, col=1)
    fig.update_xaxes(showgrid=False, row=2, col=1)
    fig.update_yaxes(zeroline=False, row=2, col=1)
    fig.update_xaxes(zeroline=False, row=2, col=1)
    if x_axis == "minmax_neg_func_score":
        fig.update_layout(title_font_family=font_list[idx_font],
                          title_text="Pathogenicity classifiers compared",
                          title_y=0.925,
                          title_font_color=yel,
                          title_font_size=18)

    return fig


# Callback allows components to interact--------------------------------------------------------------------------------------------------
@app.callback(
    Output(component_id=overview_graph, component_property='figure'),
    # Output(overview_title, 'children'),
    Input(overview_dropdown, 'value'),
    Input(overview_display, 'value'),
    Input(color_blind_option, 'on')
)
def update_overview_graph(column_name, y_axis_nucleotide,color_blind):
    df_temp = df
    c_green_red = [mid_green, '#9bbf85', '#bbd4a6', '#d7dc99', '#fff8c3', '#f8d192', '#f3a66e', '#ec785c', mid_red]
    c_blind_friendly = [mid_green, '#9bbf85', '#bbd4a6', '#e7f7d5', '#fcf2f8', '#f6d3e8', '#d091bb', '#b3589a',
                        mid_purple]
    if color_blind:
        colors = c_blind_friendly
    else:
        colors = c_green_red
    dict_color_consq = {"Neutral": colors[0], "Intermediate": colors[4], "Loss of Function": colors[-1],
                        "Synonymous": colors[0], "Intron": colors[1], "5' UTR": colors[2],
                        "Splice region": colors[3], "Missense": colors[4], "Splice acceptor": colors[5],
                        "Splice donor": colors[6], "Start lost": colors[7], "Stop gained": colors[8],
                        "Benign": colors[0], "Likely benign": colors[1],
                        "Uncertain significance": colors[3], "Absent": colors[4],
                        "Conflicting interpretations of pathogenicity": colors[5],
                        "Likely pathogenic": colors[6], "Pathogenic/Likely pathogenic": colors[7],
                        "Pathogenic": colors[8]}

    list_temp = [x for x in list(dict_color_consq.keys()) if x in list(df[column_name])]

    # Create a categorical variable with the desired order
    df_temp[column_name] = pd.Categorical(df[column_name], categories=list_temp, ordered=True)
    # Sort the dataframe based on the categorical variable
    df_temp = df_temp.sort_values(column_name)

    size_marker, height_grph, marker_symb = 8, 340, "square"
    limit = (0.5, 1.25)
    y_axis = [limit[0] + (limit[1] - limit[0]) / 2 for x in df_temp['Genomic position']]
    yaxis_dict = dict(showgrid=False, visible=False)

    if y_axis_nucleotide == "Variants expanded by nucleotide type":
        marker_symb, size_marker, y_axis, height_grph = "square", 8, 'alt_pos', 340
        yaxis_dict = dict(showgrid=False, zeroline=False, title='Nucleotide', tickvals=[0.5, 0.75, 1, 1.25],
                          ticktext=['T', 'G', 'C', 'A'])
    fig = px.scatter(data_frame=df_temp,
                     x=df_temp['Genomic position'],
                     y=y_axis,
                     height=height_grph,
                     color_discrete_map=dict_color_consq,
                     # size="1/AC",
                     color=column_name,
                     custom_data=["SGE",
                                  "UKB", 'Consequence', 'Clinvar',
                                  'cohort_allele_count', 'var_name', 'Exon'],
                     category_orders={'label': list(dict_color_consq.keys())}
                     )

    fig.update_traces(marker=dict(size=size_marker, symbol=marker_symb), selector=dict(mode="markers"),
                      hovertemplate="<br>".join([
                          "<b>%{customdata[5]}</b>",
                          "SGE: %{customdata[0]}",
                          "UKB: %{customdata[1]}",
                          "Clinvar classication: %{customdata[3]}",
                          "Consequence: %{customdata[2]}",
                          "Number of allele count in UKB: %{customdata[4]}",
                      ]))

    # add reference variant
    if y_axis_nucleotide == "Variants expanded by nucleotide type":
        ref_plot = px.scatter(data_frame=df_temp,
                              x=df_temp['Genomic position'],
                              y=df_temp['ref_pos'],
                              height=height_grph,
                              custom_data=["Genomic position", "ref", "Exon"],
                              category_orders={'label': list(dict_color_consq.keys())}
                              )
        ref_plot.update_traces(marker=dict(size=size_marker, symbol="square-open", color=yellow),
                               selector=dict(mode="markers"),
                               name='ref', hovertemplate="<br>".join([
                "<b>Reference allele</b>",
                "Position: %{customdata[0]}",
                "Nucleotide: %{customdata[1]}"
            ]))
        fig.add_trace(ref_plot.data[0])

    # add exon
    for i in range(len(exon_list)):
        start, end = exon_list[i][1] - 0.5, exon_list[i][0] + 0.5
        enlarge = 0.4
        texex = "Exon" + str(i + 1)
        pos_xanchor = 'center'
        if i + 1 in [1, 2, 3, 4, 5, 6, 7, 8, 9]:
            texex = texex + "   "
        fig.add_shape(type="rect",
                      x0=start, y0=limit[0] - enlarge, x1=end, y1=limit[1] + enlarge,
                      line=dict(color=yel_exon, width=2),
                      fillcolor=yel_exon, layer='below')
        if i + 1 in [23, 10]:
            start = start + 250
        if i + 1 in [17]:
            start = start + 150
        fig.add_annotation(x=start, y=limit[1] + 0.1,
                           text=texex,
                           showarrow=False,
                           font=dict(color=yel, family=font_list[idx_font], size=10), textangle=270,
                           xanchor=pos_xanchor, yanchor='bottom', bgcolor=dark_gray, opacity=1)
    fig.add_shape(type="rect",
                  x0=start, y0=limit[0] - enlarge - 0.2, x1=end, y1=limit[1] + enlarge,
                  line=dict(color=transparent, width=2),
                  fillcolor=transparent, layer='below')
    # add intron
    for i in range(len(intron_list)):
        start, end = intron_list[i][0] - 0.5, intron_list[i][1] + 0.5
        fig.add_shape(type='line', x0=start, x1=end,
                      y0=limit[0] + (limit[1] - limit[0]) / 2, y1=limit[0] + (limit[1] - limit[0]) / 2,
                      line=dict(color=yel_exon, width=2), layer='below')

    fig.update_layout(
        plot_bgcolor=transparent,
        paper_bgcolor=transparent,
        xaxis=dict(showgrid=False, visible=True, linecolor=None, linewidth=1, title="Genomic position"),
        yaxis=yaxis_dict,
        font_family=font_list[idx_font],
        legend=dict(orientation='h', yanchor='top', y=1.2, xanchor='left', x=0,
                    title="<b>" + column_name + ' annotation<b>', font_family=font_list[idx_font]),
        font_color=yel,
        modebar=dict(
            # "bgcolor": yellow,'activecolor':dark_gray, 'color':yel, 'bordercolor':'blue',
            bgcolor=transparent,
            activecolor=yel,
            color=yellow)
    )

    return fig.update_layout(uirevision=True)


@app.callback(
    Output('default-molecule3d-output', 'children'),
    Output(mol_viewer_colorbar, 'figure'),
    Input('dashbio-default-molecule3d', 'selectedAtomIds')
)
def show_selected_atoms(atom_ids):
    # Get a color bar
    fig = px.scatter(df, x='rna.score.1', y='rna.score.2', color='minmax_neg_func_score',
                     color_continuous_scale=['#FFFFFF', mid_red],
                     title="Comparison of RNA replicates")

    fig.update_traces(opacity=1)
    fig.update_layout(
        plot_bgcolor=dark_gray,
        paper_bgcolor=dark_gray_transp,
        xaxis=dict(showgrid=False, visible=True, linecolor=None, linewidth=2),
        yaxis=dict(showgrid=False, visible=True, linecolor=None, linewidth=2),
        height=300,
        title_font_size=18,
        title_x=0.95,
        font_color=yel,
        font_family=font_list[idx_font],
        xaxis_title="RNA score",
        yaxis_title="RNA score",
        coloraxis_colorbar=dict(
            title="SGE fct score",
            lenmode="pixels", len=200,
            yanchor="top", y=1.5, x=-0.5,
            tickvals=[0, 1],
            tickmode='array',
            ticks="outside",
            ticktext=["Neutral", "LoF"],
            ticklabelposition='outside right'
        ))

    if atom_ids is None or len(atom_ids) == 0:
        phr1 = 'Click somewhere on the protein structure to select an amino acid.'
        return phr1, fig
    for atm in atom_ids:
        print('Residue name: {}'.format(data['atoms'][atm]['residue_name']))

    aa_name = 'Reference amino acid: ', data['atoms'][atm]['residue_name'], ', position: ',\
              str(data['atoms'][atm]['residue_index']), ', chain: ', str(data['atoms'][atm]['chain']),'\n'
    subset_df = df.loc[df['aa_pos'] == data['atoms'][atm]['residue_index']]
    print(subset_df)
    if len(list(subset_df['var_name'])) < 1:
        return html.Div(
            [html.Br(), html.Div(aa_name),
             html.Div("No variant tested in SGE and present in UKB correspond to this amino acid"), html.Br()]), fig

    else:
        if len(list(subset_df['var_name'])) == 1:
            variant_nb = str(len(list(subset_df['var_name']))) + ' variant at this position.' + '\n'
        if len(list(subset_df['var_name'])) > 1:
            variant_nb = str(len(list(subset_df['var_name']))) + ' variants at this position.' + '\n'
        variant_list = []
        for i in range(len(list(subset_df['var_name']))):
            variant_list.append(str(list(subset_df['var_name'])[i]) + ', '+' SGE function score: ' + str(
                round(list(subset_df['minmax_neg_func_score'])[i],3)) + ', variant amino acid: '+ str(list(subset_df['aa_alt'])[i]))
        return html.Div([html.Br(), html.Div(aa_name),
                         html.Div(variant_nb)] +
                        [html.Div(var) for var in variant_list] +
                        [html.Br()]), fig


@app.callback(
    Output(component_id=three_d_graph, component_property='figure'),
    Output(component_id=clinvar_hist_graph_sge, component_property='figure'),
    Output(component_id=clinvar_hist_graph_ac, component_property='figure'),
    Output(component_id=clinvar_hist_graph_cadd, component_property='figure'),
    Input(component_id=overview_graph, component_property="selectedData"),
    Input(color_blind_option, 'on'),
    Input(exon_option, 'on')
)
def update_3d_graph(slct_data, color_blind, exon_option):
    fig3 = histogram("minmax_neg_func_score", color_blind)
    fig4 = histogram("cadd_score", color_blind)
    fig5 = histogram("1/AC", color_blind)
    black3dbg = dict(showbackground=True, backgroundcolor=transparent, gridcolor=light_gray, gridwidth=0.5,
                     zeroline=False)
    subtittle = "<br><sup>Choose the rectangle tool in the menu bar of the gene overview above to subset variants of interest.</sup>"
    if exon_option:
        color_exons = 'Exon'
        legend_showing = True
    else:
        color_exons = 'cumulative_score'
        legend_showing = False
    if slct_data is None or slct_data == {'points': []}:
        fig2 = px.scatter_3d(df, x='cadd_score', y='minmax_neg_func_score', z='1/AC',
                             color=color_exons, color_discrete_sequence=exons_color_l1,
                             custom_data=["var_name", 'Exon', 'cohort_allele_count', 'cadd_score',
                                          'minmax_neg_func_score'])  # color_continuous_scale=[('rgb(246, 190, 0)'),(),()])

        fig2.update_traces(hovertemplate="<br>".join([
            "<b>%{customdata[0]}</b>",
            "Exon: %{customdata[1]}",
            "SGE function score: %{customdata[4]}",
            "CADD score: %{customdata[3]}",
            "Number of allele count in UKB: %{customdata[2]}",
        ]))
        if not exon_option:
            fig2.update_traces(marker=dict(size=4, autocolorscale=True, color=df['cumulative_score']))
        else:
            fig2.update_traces(marker=dict(size=4))
        fig2.update_layout(scene=dict(
            xaxis_title=dict(text='CADD score', font=dict(color=yel)),
            yaxis_title=dict(text='SGE fct score', font=dict(color=yel)),
            zaxis_title=dict(text='1/AC', font=dict(color='rgba(248,255,24)')),
            xaxis=black3dbg,
            yaxis=black3dbg,
            zaxis=black3dbg,
            bgcolor=dark_gray),
            paper_bgcolor='rgb(41,41,41)',
            font_family=font_list[idx_font],
            font_color=yel,
            title_font_family=font_list[idx_font],
            title_text="Allele frequency, CADD and SGE function score for all variants" + subtittle,
            title_y=0.955,
            title_font_color=yel,
            title_font_size=18,
            legend=dict(orientation='v', yanchor='top', y=0.9, xanchor='left', x=0, title="Region",
                        font=dict(color=yel)),
            showlegend=legend_showing,
            height=813, coloraxis_colorbar=dict(
                title="Cumulative score",
                lenmode="pixels", len=380, thickness=8
            ))

        return fig2, fig3, fig4, fig5

    if slct_data['points'] == []:
        fig2 = px.scatter_3d(title="Please select at least one variant")
        fig2.update_layout(scene=dict(
            bgcolor=dark_gray),
            paper_bgcolor=dark_gray,
            font_family=font_list[idx_font],
            font_color=yel,
            title_font_family=font_list[idx_font],
            title_font_color=yel,
            title_font_size=18,
            height=785)
        return fig2, fig3, fig4, fig5
    else:
        print(f'selected data: {slct_data}')
        exons = [slct_data['points'][i]['customdata'][-1] for i in range(len(slct_data['points']))]
        var = [slct_data['points'][i]['customdata'][-2] for i in range(len(slct_data['points']))]

        dff2 = df[df.var_name.isin(var)]
        fig2 = px.scatter_3d(dff2, x='cadd_score', y='minmax_neg_func_score', z='1/AC',
                             color=color_exons,
                             custom_data=["var_name", 'Exon', 'cohort_allele_count', 'cadd_score',
                                          'minmax_neg_func_score'],
                             color_discrete_sequence=exons_color_l1)  # \nin "+str(set(exons)).replace("{", '').replace("'", '').replace("}", ''))

        # add subptitkes

        fig2.update_traces(hovertemplate="<br>".join([
            "<b>%{customdata[0]}</b>",
            "Exon: %{customdata[1]}",
            "SGE function score: %{customdata[4]}",
            "CADD score: %{customdata[3]}",
            "Number of allele count in UKB: %{customdata[2]}",
        ]))
        if not exon_option:
            fig2.update_traces(marker=dict(size=4, autocolorscale=True, color=df['cumulative_score']))
        else:
            fig2.update_traces(marker=dict(size=4))

        fig2.update_layout(scene=dict(
            xaxis_title=dict(text='CADD score', font=dict(color=yel)),
            yaxis_title=dict(text='SGE fct score', font=dict(color=yel)),
            zaxis_title=dict(text='1/AC', font=dict(color='rgba(248,255,24)')),
            xaxis=black3dbg,
            yaxis=black3dbg,
            zaxis=black3dbg,
            bgcolor=dark_gray),
            paper_bgcolor='rgb(41,41,41)',
            font_family=font_list[idx_font],
            font_color=yel,
            title_font_family=font_list[idx_font],
            title_text="Allele frequency, CADD and SGE function score for all variants",
            title_y=0.955,
            title_font_color=yel,
            title_font_size=18,
            legend=dict(orientation='v', yanchor='top', y=0.9, xanchor='left', x=0, title="Region",
                        font=dict(color=yel)),
            showlegend=legend_showing,
            height=785, coloraxis_colorbar=dict(
                title="Cumulative score",
                lenmode="pixels", len=380, thickness=8
            ))

        return fig2.update_layout(uirevision=True), fig3, fig4, fig5


# Run app
if __name__ == '__main__':
    app.run_server(debug=True)