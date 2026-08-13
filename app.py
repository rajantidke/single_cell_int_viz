#!/usr/bin/env python3

import dash
import dash_bootstrap_components as dbc
from dash import html, page_container

external_stylesheets = [
    dbc.themes.DARKLY,
]

app = dash.Dash(
    __name__,
    use_pages=True,
    external_stylesheets=external_stylesheets,
    suppress_callback_exceptions=True,
)
server = app.server

app.layout = html.Div([
    page_container
])

if __name__ == '__main__':
    app.run(debug=False)