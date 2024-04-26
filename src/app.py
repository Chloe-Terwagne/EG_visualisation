"""
    Project :
    Description:
    Name : Chloé Terwagne
    date : 24 April 2024
    Python version 3.10
"""
# IMPORT ---------------------------------------------------------------

import numpy as np
from dash import Dash, dcc, html, Output, Input
import dash_bootstrap_components as dbc
import plotly.express as px
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots


# FONT & COLOR  ---------------------------------------------------------------
font_list = ["Arial", "Balto", "Courier New", "Droid Sans", "Droid Serif", "Droid Sans Mono", "Gravitas One",
             "Old Standard TT", "Open Sans", "Overpass", "PT Sans Narrow", "Raleway", "Times New Roman"]
idx_font = -1  # UCL font
# Color Model
pink_vibrant = 'rgb(172, 20, 90)'  # UCL color
blue_vibrant = 'rgb(52, 198, 198)'  # UCL color
yellow_vibrant = 'rgb(255, 202, 56,0.5)'  # UCL color

yellow_mutated = 'rgb(218, 214, 202)'  # UCL color
pink_mutated = 'rgb(222, 184, 195)'  # UCL color
purple_mutated = "rgb(198,176,188)"  # UCL color
transparent = 'rgba(0,0,0,0)'


# pd display
pd.set_option('display.width', 900)
pd.set_option('display.max_columns', 200)
pd.set_option("display.max_rows", None)


# FUNCTION  -----------------------------------------------------------------------------------------------------------
def preprocess_df(df):
    #df.dropna(subset=['penetrance'], inplace=True)
    # groupd disease with fw counts into a others column
    group_counts = df['diseaseGroup'].value_counts()
    # Identify disease groups with less than 30 counts
    groups_to_replace = group_counts[group_counts < 85].index
    # Replace these groups with "other"
    df.loc[df['diseaseGroup'].isin(groups_to_replace), 'diseaseGroup'] = 'Other'
    df = df[df['diseaseGroup'] != 'Other']

    # Count the occurrences of each disease group
    group_counts = df['diseaseGroup'].value_counts()

    # Sort the disease groups by decreasing count
    sorted_groups = group_counts.sort_values(ascending=False)
    df['penetrance'] = df['penetrance'].replace(['nan', 'unknown'], 'Unknown')
    # Dictionary to map old column names to new column names
    column_name_mapping = {"strand_y":'DNA strand',
     "penetrance": 'Penetrance',
     "CEGs": 'Core-essential genes',
     }
    # Rename columns using the rename() method
    df.rename(columns=column_name_mapping, inplace=True)
    df["Confidence level"]=df["confidenceLevel"].map({1.0:'Low',2.0:'Intermediate',3.0:'High' })

    return df, sorted_groups

def get_chr_dict(df):
    # Assuming df is your DataFrame containing the relevant data
    # Subset the DataFrame to unique 'geneSymbol_panel_app_x' column values
    df_unique_genes = df[['gene_name', 'seqname']].drop_duplicates()
    df_unique_genes = df_unique_genes.dropna()

    # Initialize a dictionary to store the count of rows by 'seqname' for each unique gene
    dct_chr = {}
    for chr in df_unique_genes["seqname"].unique():
        if chr != np.nan:
            dct_chr[chr]=df_unique_genes[df_unique_genes["seqname"]==chr].shape[0]
    return dct_chr

# MAIN ---------------------------------------------------------------------------------------------------------------

# cleaning data
text_b1 = "This study is a crucial component of my analysis aimed at prioritizing a large number of variants for editing via CRISPR-Cas9 techniques. The final selection encompasses mutations across multiple genes, making the identification of suitable genes for the cell line I'm working with an essential and critical step."
text_b2 ="Only genes genes essential for HAP1 survival and found in at least one gene panel are displayed (N=631). It's noteworthy that genes may appear multiple times for each condition, as each dot represents a unique gene-condition pair. The HAP1 essentiality ratio utilized in this analysis is sourced from Blomen et al. (2015). Gene-disease associations are established through the Gene Panel App of Genomics England. Conserved regions are retrieved from gnomAD v4, while pathogenic and likely pathogenic submissions to ClinVar are gathered as of 24/04/2024."

df = pd.read_csv("https://github.com/Chloe-Terwagne/EG_visualisation/blob/main/data/output_df/20240425_5226_genes_disease_selected.tsv?raw=true", sep='\t')

df, sorted_groups = preprocess_df(df)
chromosome_data=get_chr_dict(df)

# Launch app------------------------------------------------------------------------------------------------
app = Dash(__name__, suppress_callback_exceptions=True,external_stylesheets=[dbc.themes.BOOTSTRAP],
           meta_tags=[{'name': 'viewport', 'content': 'width=device-width, initial-scale=1.0'}])
server = app.server

# Build your components------------------------------------------------------------------------------------------------

overview_title = dcc.Markdown(children='', style=dict(font_family=font_list[idx_font], font_color=pink_vibrant))
overview_display = dcc.RadioItems([
    {
        "label": html.Div(['Core-essential genes'], style={'color':pink_vibrant, 'font-size': 15}),
        "value": "Core-essential genes",
    },
    {
        "label": html.Div(['DNA strand'], style={'color': pink_vibrant, 'font-size': 15}),
        "value": "DNA strand",
    },
    {
        "label": html.Div(['Confidence level'], style={'color': pink_vibrant, 'font-size': 15}),
        "value": "Confidence level",
    },
    {
        "label": html.Div(['Penetrance'], style={'color': pink_vibrant, 'font-size': 15}),
        "value": "Penetrance",
    },
], value='Confidence level', inline=True, labelStyle={"margin-right": "30px"})

pie_graph = dcc.Graph(figure={},  config={
        'staticPlot': False,
        'scrollZoom': False,
        'doubleClick': 'reset',
        'showTips': True,
        'displayModeBar': False,  # Set to False to hide the ModeBar
        'displaylogo': False,
        'modeBarButtonsToRemove': ['lasso2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d'],
    }, selectedData=None)
overview_graph = dcc.Graph(figure={}, config={'staticPlot': False, 'scrollZoom': False, 'doubleClick': 'reset',
                                              'showTips': True,'displayModeBar': 'hover', 'displaylogo': False,
                                              'modeBarButtonsToRemove': ['lasso2d', 'zoomIn2d', 'zoomOut2d',
                                                                         'autoScale2d'], }, selectedData=None)
# Custom sorting function to order keys
ordered_chromosome_data={'1': 69,'2': 44, '3': 45, '4': 23, '5': 25, '6': 24, '7': 23, '8': 15, '9': 24, '10': 19, '11': 30, '12': 32, '13': 9, '14': 23, '15': 13, '16': 40, '17': 51, '18': 12, '19': 35,  '20': 16, '21': 9, '22': 12, 'X': 36}
label= list(ordered_chromosome_data.keys())
bar_plot = dcc.Graph(
    figure={
        'data': [
            go.Bar(
                x=list(range(0,24)),
                y=list(ordered_chromosome_data.values()),
                marker_color=yellow_mutated,
                hoverinfo='y+text',  # Enable hovering to show y value
                # text=[f'Count: {count}' for count in ordered_chromosome_data.values()],  # Text to show on hover

            )
        ],
        'layout': {
            'title': 'Number of essential genes by chromosome',
            'plot_bgcolor' : transparent,
            'paper_bgcolor' : transparent,
            'xaxis': {
                'title': '',
                'showticklabels': True,  # Show tick labels on x-axis
                'tickvals': list(range(0,24)),  # Set tick positions
                'ticktext': list(ordered_chromosome_data.keys()),  # Set tick labels
                'tickangle': 0,  # Set tick label angle to horizontal
                'font_size':2,

            },
            'font':dict(family=font_list[idx_font], color=pink_vibrant),
            'font_color' : yellow_mutated,
            'height': 300,
            'yaxis': {'visible': False},  # Remove y-axis
            'margin': {'l': 8, 'r': 10, 't': 25, 'b': 25}  # Set margin to zero
        },
    },    config={'displayModeBar': False}  # Set to False to hide the ModeBar

)
# Customize Layout--------------------------------------------------------------------------------------------

app.layout = \
    dbc.Container([
        dbc.Row([],  style={ 'padding': '40px', 'margin-top':'0px', 'backgroundColor': yellow_mutated}),
        dbc.Row([
            dbc.Col([
                dbc.Card(
                    [html.Div(children=[
                                html.H1("Unveiling Gene Essentiality for Disease Insights", style={'color': pink_vibrant,
                                                                                      'font-size':28,
                                                                                        'text-align': 'center',
                                                                                        'margin-top': '30px',
                                                                                                   'margin-right': '30px',

                                                                                                   'margin-bottom': '20px',
                                                                                        'font-weight': 'bold',
                                                                                        'font-family': font_list[0]}),
                                html.P('In this study, we delve into genes crucial for the Haploid human cell line, HAP1, aiming to uncover their mechanisms and disease implications. Through genetic manipulations such as knockouts or variant insertions, we anticipate observing pronounced phenotypic effects in the lab setting. Exploring these complex data helps to prioritize genes that could be efficiently tested in the lab, facilitating targeted investigations into their functional roles and potential therapeutic relevance.'),
                                html.P('The plot on the right showcases genes deemed vital for HAP1 survival, where a lower ratio underscores greater gene essentiality. These genes are grouped into the top eight disease categories, ranked from the group with the highest gene count to the lowest.'),

                                # Here add a static barplot
                                html.Br(),
                                bar_plot,
                                html.Br(),
                                html.P(
                                    'The gene-disease group associations are derived from a meticulously curated effort by experts in each condition, reflecting various confidence levels denoted by plot color. Furthermore, the plot can also be color-coded to highlight core-essential genes, which are crucial in 80% or more tested cell lines according to the OGEE V3 database, as well as by the DNA strand where the gene is located and disease penetrance. Additionally, the size of each dot corresponds to the number of pathogenic or likely pathogenic variants reported in the ClinVar clinical database. The pie chart above each group illustrates their rank and the diversity of sub-disease groups within.'),
                            ])
                    ], style={'overflow-y': 'auto', 'max-height':'1200px', 'background-color': purple_mutated, 'padding': '20px','border': 'none'}
                )], style={'background-color': purple_mutated, 'padding': '20px', 'margin-left':'30px'}, width=4),
            dbc.Col([
                    dbc.Row([pie_graph],style={'margin': '0px','padding-left': '67px','padding-right': '22px', 'z-index': '1'} ),
                    dbc.Row([overview_graph],style={'padding': '0px', 'margin': '0px' , 'z-index': '0'})]),]),
        dbc.Row([html.Br()]),
        dbc.Row([overview_display],style={'padding-bottom': '80px', 'margin-left': '10px'}),
        dbc.Row([html.P(text_b2 + text_b1,style= {'color': 'rgb(128,128,128)',
                                                        'font-size':14,
                                                        'text-align': 'center'})], style={'padding': '50px', 'margin': '0px', 'backgroundColor': yellow_mutated, }),

    ], fluid=True, style={'backgroundColor': yellow_mutated})

# Callback allows components to interact--------------------------------------------------------------------------------------------------

@app.callback(
    Output(component_id=overview_graph, component_property='figure'),
    Input(overview_display, 'value'),
)
def update_overview_graph(column_name_color):
    colors = [yellow_vibrant, blue_vibrant, pink_vibrant]
    df_temp = df.copy()
    print('sorted_groups', sorted_groups)
    # get very small value for nan -> bc probably 0
    df_temp['Alleles_reported_Pathogenic_Likely_pathogenic'].fillna(0.1, inplace=True)
    df_temp.fillna("N/A", inplace=True)



    # Create a rank column based on the count of disease groups
    rank_dict = {group: rank for rank, group in enumerate(sorted_groups.index, 1)}
    df_temp['diseaseGroup_rank'] = df_temp['diseaseGroup'].map(rank_dict)

    # Add jitter effect to the rank values
    jitter_amount = 0.08  # Adjust this value as needed
    jittered_indices = np.random.uniform(-jitter_amount, jitter_amount, len(df_temp))
    df_temp['diseaseGroup_rank_jittered'] = df_temp['diseaseGroup_rank'] + jittered_indices

    # Create the scatter plot
    fig = px.scatter(df_temp, x='diseaseGroup_rank_jittered', y='ratio',
                     labels={'diseaseGroup_rank_jittered': 'Disease Group Rank', 'ratio': 'Ratio'},
                     color=column_name_color,
                     color_discrete_sequence=colors,
                     custom_data=["entityName",'ENSEMBL_ID','Confidence level' ,"omimGene", 'Alleles_reported_Pathogenic_Likely_pathogenic', 'publications', 'modeOfInheritance', 'mis.oe']
                     )

    # Update x-axis tick labels
    fig.update_xaxes(tickvals=df_temp['diseaseGroup_rank'], ticktext=df_temp['diseaseGroup'])

    # Define scaling factors
    linear_threshold = 100  # Adjust this threshold as needed
    linear_factor = 0.25  # Adjust this factor as needed
    logarithmic_factor = 0.05  # Adjust this factor as needed

    # Scale the marker size based on the column values
    marker_size_linear = np.where(df_temp['Alleles_reported_Pathogenic_Likely_pathogenic'] <= linear_threshold,
                           df_temp['Alleles_reported_Pathogenic_Likely_pathogenic'] * linear_factor,
                           np.log(df_temp['Alleles_reported_Pathogenic_Likely_pathogenic'] + 1) * logarithmic_factor)
    min_marker_size = 8.5  # Adjust this value as needed
    # Apply the minimum size constraint
    marker_size = np.maximum(marker_size_linear, min_marker_size)
    # Update traces with scaled marker size
    fig.update_traces(marker=dict(size=marker_size, symbol="circle"), selector=dict(mode="markers"),
                      hovertemplate="<br>".join([
                          "<b>%{customdata[0]}</b>",
                          "Ensembl ID: %{customdata[1]}", "Disease association confidence: %{customdata[2]}",
                          "omim Gene: %{customdata[3]}", "Publications: %{customdata[5]}",
                          "Mode of inheritance: %{customdata[6]}", "o/e missense mutation rate: %{customdata[7]}",
                          "Number Alleles reported Pathogenic: %{customdata[4]}"
                      ]))

    # Set x-axis tick values and labels
    tickvals = df_temp['diseaseGroup_rank']
    ticktext = df_temp['diseaseGroup'].str.replace('Programme', '')  # Replace spaces with line breaks
    ticktext = ticktext.str.replace('Dysmorphic and congenital', 'Congenital')  # Replace spaces with line breaks
    ticktext = ticktext.str.replace('Dysmorphic and congenital', 'Congenital')  # Replace spaces with line breaks
    ticktext = ticktext.str.replace('Neurology and', 'Neurology_&')  # Replace spaces with line breaks

    ticktext = ticktext.str.replace(' ', '<br>')  # Replace spaces with line breaks
    ticktext = ticktext.str.replace('_', ' ')  # Replace spaces with line breaks

    fig.update_layout(
        margin=dict(t=0, b=-0, l=0, r=0),  # Set all margins to 0
        plot_bgcolor=transparent,
        paper_bgcolor=transparent,
        yaxis=dict(showgrid=False, visible=True, linecolor=transparent, linewidth=1, title="HAP1 Gene Essentiality", range=[0,0.46]), font_size=18,
        font_family=font_list[idx_font],
        legend=dict(
            orientation='h',  # Set orientation to horizontal
            yanchor='bottom',  # Anchor the legend to the bottom of the plot
            y=-0.3,  # Position the legend slightly above the bottom
            xanchor='left',  # Anchor the legend to the left of the plot
            #x=0,  # Position the legend at the leftmost position
            title="<b>" + column_name_color + ' annotation<b>',  # Set legend title
            font_family=font_list[idx_font]  # Set font family for legend
        ),
        xaxis=dict(
            tickmode='array',
            tickvals=tickvals,
            ticktext=ticktext,
            showgrid=False,
            visible=True,
            linecolor=transparent,
            linewidth=1,
            title=""
        ),
        font_color=pink_vibrant,
        modebar=dict(
            # "bgcolor": yellow,'activecolor':dark_gray, 'color':yel, 'bordercolor':'blue',
            bgcolor=transparent,
            activecolor=blue_vibrant,
            color=blue_vibrant),
        height=900,  #
    )

    return fig.update_layout(uirevision=True)


@app.callback(
    Output(component_id=pie_graph, component_property='figure'),
    Input(overview_display, 'value'),
)
def update_pie(column_name_color):
    df_temp = df.copy()

    sorted_group_idx= sorted_groups.index.tolist()
    # Create subplots, using 'domain' type for pie charts
    num_groups = len(sorted_group_idx)
    specs = [[{'type': 'domain'}] * num_groups]  # 1 row, num_groups columns
    fig = make_subplots(rows=1, cols=num_groups, specs=specs, vertical_spacing = 0)
    corlors = ['rgb(128, 0, 0)', 'rgb(165, 42, 42)', 'rgb(160, 82, 45)', 'rgb(139, 69, 19)', 'rgb(210, 105, 30)',
               'rgb(205, 133, 63)', 'rgb(184, 134, 11)', 'rgb(218, 165, 32)', 'rgb(244, 164, 96)',
               'rgb(188, 143, 143)', 'rgb(210, 180, 180)','rgb(222, 184, 135)']
    corlors = [pink_mutated, purple_mutated, 'rgb(201,209,168)']
    # Define pie charts
    for i, group in enumerate(sorted_group_idx):
        # Filter data for the current disease group
        group_data = df_temp[df_temp['diseaseGroup'] == group]

        # Create a pie plot for subgroup proportions
        pie_fig = go.Pie(labels=group_data['diseaseSubgroup'],
                         values=group_data['ratio'],
                         name=group,
                         marker=dict(colors=corlors))

        # Calculate the center position of the pie
        start_pos = -0.0232
        end_pos = 1.025
        pie_center = start_pos + (end_pos - start_pos) * (2 * i + 1) / (2 * num_groups)

        # Add pie plot to the subplot with the calculated domain
        fig.add_trace(pie_fig, row=1, col=i + 1)

        # Add annotation at the center of the pie
        fig.add_annotation(
            text="<b>"+str(i + 1)+"</b>", x=pie_center, y=0.5, font_size=30, showarrow=False,
            font=dict(family=font_list[idx_font], color=pink_vibrant)
        )

    # Tune layout and hover info
    fig.update_traces(hole=.55, hoverinfo='label+percent', textinfo='none')
    fig.update(layout_showlegend=False)

    fig.update_layout(
        plot_bgcolor=transparent,
        paper_bgcolor=transparent,
        colorway=corlors,
        height=180,  # Adjust the height as needed
        title=None,
        margin=dict(t=0, b=-0, l=0, r=0))

    # fig = go.Figure(fig)
    return fig.update_layout(uirevision=True)


# Run app
if __name__ == '__main__':
    app.run_server(debug=True, port=8051)
