# -*- coding: utf-8 -*-

# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.express as px
import pandas as pd

import numpy as np

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

df = pd.DataFrame(np.array(np.transpose([[-4,1,2],[-1,4,-5]])),
                  index=['DYCOMS','BOMEX','LBA'],
                  columns= ['defBias', 'defBiasApprox'])

#print(df)

fig = px.line(df, x=df.index, y=df.columns, labels={"df.columns":"Bias"},
              title = 'Biases of default simulation and approximations thereof.')

fig.update_yaxes(title="Bias")
fig.update_xaxes(title="Metric and region")

app.layout = html.Div(children=[
        html.H1(children='Hello Dash'),

        html.Div(children='''
                         Dash: A web application framework for Python.
                     '''),

        dcc.Graph(
                    id='example-graph',
                    figure=fig
                )
])

if __name__ == '__main__':
        app.run_server(debug=True)
