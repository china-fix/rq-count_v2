import base64
import datetime
import io

import dash
from dash.dependencies import Input, Output, State
from dash import dcc
from dash import html
from dash import dash_table
import plotly.express as px

import pandas as pd


stylesheets = ['bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=stylesheets,
                suppress_callback_exceptions=True)

app.layout = html.Div([ 
    dcc.Upload(
        id='upload-data',
        children=html.Div([
            'Drag and Drop or ',
            html.A('Select Files')
        ]),
        style={
            'width': '100%',
            'height': '60px',
            'lineHeight': '60px',
            'borderWidth': '1px',
            'borderStyle': 'dashed',
            'borderRadius': '5px',
            'textAlign': 'center',
            'margin': '10px'
        },
        # Allow multiple files to be uploaded
        multiple=True
    ),
    html.Div(id='output-div'),
    html.Div(id='output-datatable'),
])


def parse_contents(contents, filename, date):
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    try:
        if 'tab' in filename:
            # Assume that the user uploaded a tab file
            df = pd.read_csv(
                io.StringIO(decoded.decode('utf-8')), sep='\t')
    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])

    return html.Div([
        html.H5(filename),
        html.H6(datetime.datetime.fromtimestamp(date)),
        html.P("hi Xiao please Inset X axis data"),
        dcc.Dropdown(id='xaxis-data',
                     options=[{'label':x, 'value':x} for x in df.columns]),
        html.P("hi Xiao please Inset Y axis data"),
        dcc.Dropdown(id='yaxis-data',
                     options=[{'label':x, 'value':x} for x in df.columns]),
        html.Button(id="submit-button", children="Create Graph"),
        html.P(),
        dcc.Input(id="nfrom", type="number", placeholder="1-index caculate start"),
        dcc.Input(id="nend", type="number", placeholder="1-index caculate end"),
        html.Button(id="caculate-button", children="caculate"),
        html.Hr(),
        html.Div(id="number-out"),
        html.Hr(),
        dash_table.DataTable(
            data=df.to_dict('records'),
            columns=[{'name': i, 'id': i} for i in df.columns],
            page_size=15
        ),
        dcc.Store(id='stored-data', data=df.to_dict('records')),
        html.Hr(),  # horizontal line

        # For debugging, display the raw contents provided by the web browser
        html.Div('Raw Content'),
        html.Pre(contents[0:200] + '...', style={
            'whiteSpace': 'pre-wrap',
            'wordBreak': 'break-all'
        })
    ])


@app.callback(Output('output-datatable', 'children'),
              Input('upload-data', 'contents'),
              State('upload-data', 'filename'),
              State('upload-data', 'last_modified'))
def update_output(list_of_contents, list_of_names, list_of_dates):
    if list_of_contents is not None:
        children = [
            parse_contents(c, n, d) for c, n, d in
            zip(list_of_contents, list_of_names, list_of_dates)]
        return children


@app.callback(Output('output-div', 'children'),
              Input('submit-button','n_clicks'),
              State('stored-data','data'),
              State('xaxis-data','value'),
              State('yaxis-data', 'value'))
def make_graphs(n, data, x_data, y_data):
    if n is None:
        return dash.no_update
    else:
        bar_fig = px.bar(data, x=x_data, y=y_data)
        # print(data)
        return dcc.Graph(figure=bar_fig)

    
@app.callback(
    Output("number-out", "children"),
    Input('caculate-button','n_clicks'),
    Input("nfrom", "value"),
    Input("nend", "value"),
    State('stored-data','data'))
def mean_caculate(n, fval, tval, data):
    if n is None:
        return dash.no_update
    else:
        df=pd.DataFrame.from_dict(data)
        mean=df['coverage'][fval-1:tval-1].mean()
        std=df['coverage'][fval-1:tval-1].std()
        return "caculate from : {}bp to {}bp,  mean coverage={}  standard deviation={}".format(fval, tval, mean, std)
        



if __name__ == '__main__':
    app.run_server(host='0.0.0.0', port='9999', debug=False)