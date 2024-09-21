import os
import io
import math
import base64 
from dotenv import load_dotenv, find_dotenv
from PIL import Image, ImageDraw
from Bio import Entrez, SeqIO
from dash import Dash, html, dcc, callback, Input, Output, State
import dash_bootstrap_components as dbc

load_dotenv(find_dotenv())


# API Credentials and Redirect. 
client_id = os.getenv("entrez_email")
max_seq_len = int(os.getenv("max_seq_len")) #2.5M 
covid_accession = 'NC_045512'

color_dict = {"A": "Green", "T": "Red", "C": "Orange", "G": "Blue"}

# Dash App
app = Dash(__name__, external_stylesheets=[dbc.themes.LUMEN])
server = app.server
app.title = "Seq Grid"

# Layout
app.layout = html.Div([
    html.H1("Seq-Grid", style={"textAlign": "center"}),
    html.Br(),
    html.P("A fun way to visualize DNA (and in the future, Protein) sequences!", style={"textAlign": "center"}),
    html.Br(),
    dbc.Row([
        dbc.Col([
            dbc.Row([
                html.Div([
                    html.Label("Sequence or Accession Number:"),
                    dbc.Textarea(id = 'seq_input',  
                                 value = 'NC_045512', 
                                 maxLength = max_seq_len, 
                                 style = {'width': '100%'}),
                ])

            ],style={'padding': '10px'}),
            dbc.Row([
                html.Div([
                    html.Label("Shape"),
                    dcc.Dropdown(id='shape_dropdown', 
                                 options=[
                                     {'label': 'Circle', 'value': 'circle'},
                                     {'label': 'Square', 'value': 'square'}
                                 ],
                                 value='circle'
                                 ),                    
                    
                ], style={'padding': '10px'})
            ],style={'padding': '10px'}),
            dbc.Row([
                dbc.Col([
                    html.Label("Shape Size:"),
                    dbc.Input(id='shape_size', type='number', value=5, min=2),
                ]),
                dbc.Col([
                    html.Label("Shape Pad:"),
                    dbc.Input(id='shape_pad', type='number', value=2, min=0),
                ]),
            ],style={'padding': '10px'}),
            dbc.Row([
                dbc.Col([
                    html.Label("X Ratio:"),
                    dbc.Input(id='x_ratio', type='number', value=1, min=1)
                    ]),
                dbc.Col([
                    html.Label("Y Ratio:"),
                    dbc.Input(id='y_ratio', type='number', value=1, min=1),
                ]),
            ],style={'padding': '10px'}),
            dbc.Row([
                dbc.Col([
                    html.Label("External Margin:"),
                    dbc.Input(id='ext_margin', type='number', value=5, min=0),
                ]),
                dbc.Col([
                    html.Label("Tight Ratio:"),
                    dbc.Checkbox(id='tight_ratio', value=False),
                ]),
            ],style={'padding': '10px'}),
            dbc.Row([
                html.Div([
                    dbc.Button("Generate Image", id='submit', n_clicks=0, color="primary"),
                ]),
                
            ],style={'padding': '10px'}),
            html.Br(),
            html.Div(id='img_dwnld_btn_div', children=[
                dbc.Button("Download Image", id='img_dwnld_btn', n_clicks=0, color="primary")
            ], style={'padding': '10px', 'display':'none'}),
        ], width=3),
        dcc.Download(id="download_holder"),
        dbc.Col([
            html.Div(id='image_display', 
                     style={'border': '1px solid #ccc', 'height': '100%', 'width': '100%'}, 
                     children=[])            
        ], width=9),
    ]),


    html.Div([
        "Made by ", 
        html.A("James Welch", href="https://github.com/jwelch1123",  target="_blank", style={'color': 'grey'}),
        " · ",
        "View ", 
        html.A("Seq-Grid on GitHub", href="https://github.com/jwelch1123/seq-grid",  target="_blank", style={'color': 'grey'})
        ], style={'textAlign': 'center', 'color': 'grey', 'fontSize': '0.8em', 'marginTop': '20px'}
        )


    ],
    style={'margin': 'auto',
           'marginTop': '1%',
           'maxWidth': '95%', 
           'maxHeight': '95%',
           'padding': '20px',
           'border':'1px solid #ccc',
           'borderRadius': '10px'
           }
)

# Functions
def find_ratio_factors(l: int, x:int, y:int, tight: bool = False):
    '''
        Returns the scale factor for x & y of the smallest grid large enough to contain n items 
    '''
    scale_factor = 1
    while (x * y * (scale_factor ** 2)) < l:
        scale_factor += 1
        if scale_factor > l:
            break

    x = x * scale_factor
    y = y * scale_factor

    if tight:
        y = math.ceil(l / x)
        x = l if l < x else x

    return x, y 

def get_genbank_sequence(accession: str):
    Entrez.email = os.getenv("entrez_email")
    handle = Entrez.efetch(
                        db="nucleotide", 
                        id=accession, 
                        rettype="gb", 
                        retmode="text")
    try:
        record = SeqIO.read(handle, "genbank")
        handle.close()
        return record.seq
    except:
        return ''
    
def draw_grid(dims: tuple, seq: list, shape: str, shape_size: int, shape_pad: int, ext_margin: int, color_dict: dict, bkgd_color: str = "black"):
    '''
        Draws the grid
    '''

    img_width  = (ext_margin * 2) + (dims[0] * shape_size) + ((dims[0] - 1) * shape_pad) 
    img_height = (ext_margin * 2) + (dims[1] * shape_size) + ((dims[1] - 1) * shape_pad)
    
    img = Image.new('RGB', (img_width, img_height), color = bkgd_color) 
    draw = ImageDraw.Draw(img)

    for i, val in enumerate(seq):

        x = i % dims[0] 
        y = i // dims[0]
        
        x_coord = ext_margin + (x * (shape_size + shape_pad))
        y_coord = ext_margin + (y * (shape_size + shape_pad))

        x2 = x_coord + shape_size-1
        y2 = y_coord + shape_size-1

        color = color_dict[val]

        #draw.rectangle([x_coord, y_coord, x2, y2], fill = color)
        #draw.circle(((x_coord+x2)/2, (y_coord+y2)/2),radius=(x2-x_coord)/2, fill = color)
        create_shape(draw, shape, x_coord, y_coord, x2, y2, color)

    return img
        
def pil_to_b64(im, enc_format="png", **kwargs):
    """
    Converts a PIL Image into base64 string for HTML displaying
    :param im: PIL Image object
    :param enc_format: The image format for displaying. If saved the image will have that extension.
    :return: base64 encoding
    """
    buff = io.BytesIO()
    im.save(buff, format=enc_format, **kwargs)
    encoded = base64.b64encode(buff.getvalue()).decode("utf-8")

    return encoded

def create_shape(draw, shape: str, x_coord: int, y_coord: int, x2: int, y2: int, color: str):
    if shape == "square":
        return draw.rectangle([x_coord, y_coord, x2, y2], fill = color)
    elif shape == "circle":
        return draw.circle(((x_coord+x2)/2, (y_coord+y2)/2),radius=(x2-x_coord)/2, fill = color)



# Callbacks
@app.callback(
        Output('image_display', 'children'), 
        Input('submit', 'n_clicks'),
        State('seq_input', 'value'),
        State('shape_dropdown', 'value'),
        State('shape_size', 'value'),
        State('shape_pad', 'value'),
        State('x_ratio', 'value'),
        State('y_ratio', 'value'),
        State('ext_margin', 'value'),
        State('tight_ratio', 'value')
        )
def create_n_display(n_clicks, seq_input, shape_dropdown, shape_size, shape_pad, x_ratio, y_ratio, ext_margin, tight_ratio):
    # Need to do validation
    if n_clicks == 0:
        return html.Img(style={'width': '100%'})
    elif not seq_input:
        return html.P("Please enter a sequence or accession number.")
    elif len(seq_input) > max_seq_len:
        return html.P(f"Sequence length exceeds maximum length of {max_seq_len}")
    

    clean_seq = ''.join(filter(str.isalpha, seq_input)).strip().upper()
    if not all([char in color_dict.keys() for char in clean_seq]):
        genbank_cache = get_genbank_sequence(seq_input)
        if genbank_cache:
            clean_seq = genbank_cache
        else:
            invalids = {[char for char in clean_seq if char not in color_dict.keys()]}
            return html.P(f"Invalid sequence or accession number. Please try again.")


    # try:
    #     clean_seq = get_genbank_sequence(seq_input)
    # except:
    #     clean_seq = ''.join(filter(str.isalpha, seq_input)).strip().upper()

    #     if not all([char in color_dict.keys() for char in clean_seq]):
    #         invalids = {[char for char in clean_seq if char not in color_dict.keys()]}
    #         return html.P(f"Invalid sequence. Please enter a valid DNA sequence. Invalid characters: {', '.join(invalids)}")
    #     else:
    #         return html.P("Invalid Accession Number. Please enter a valid accession number or a DNA sequence.")


    #clean_seq = ''.join(filter(str.isalpha, seq_input)).strip().upper()

    #if not all([char in color_dict.keys() for char in clean_seq]):
    #    invalids = {[char for char in clean_seq if char not in color_dict.keys()]}
    #    return html.P(f"Invalid sequence. Please enter a valid DNA sequence. Invalid characters: {', '.join(invalids)}")
    #else:    


    x, y = find_ratio_factors(len(clean_seq), x_ratio, y_ratio, tight_ratio)
    
    img = draw_grid((x,y), clean_seq, shape_dropdown, shape_size, shape_pad, ext_margin, color_dict)
    img_encoded = 'data:image/png;base64,' + pil_to_b64(img)
    
    return html.Img(src=img_encoded, style={'width': '100%', 'image-rendering':'pixelated'})


@app.callback(
        Output('img_dwnld_btn_div', 'style'),
        Input('image_display', 'children')
)
def create_download_btn(img):
    try:
        if img['props']['src']:
            return {'padding': '10px'}
        else:
            return {'padding': '10px', 'display':'none'}
    except:
        return {'padding': '10px', 'display':'none'}
    
@app.callback(
        Output('download_holder','data'),
        Input('img_dwnld_btn', 'n_clicks'),
        State('image_display', 'children'))
def download_image(n_clicks, image_div_children):
    if image_div_children and (n_clicks > 0):
        img_src = image_div_children['props']['src'].split(',')[1]
        img_bytes = base64.b64decode(img_src)
        return dcc.send_bytes(img_bytes, filename='seq-grid.png')
        
# Do Things
if __name__ == '__main__':
    app.run_server(debug=True)