import dash
from dash import html, dcc

dash.register_page(__name__, path='/')

# Personal info
TEAM_INFO = [
    {
        "name": "Rajan Ranjeet Tidke",
        "github": "https://github.com/rajantidke",
        "email": "mailto:rajan.tidke.7@gmail.com"
    }
]

# Function to create team member entries with icons
def create_team_member_entry(member):
    return html.Div([
        # Member name
        html.Span(member["name"], style={
            'color': '#ff8aff',
            'fontSize': '1.2vw',
            'fontFamily': "'DM Mono', monospace",
            'marginRight': '0.8vw'
        }),

        # GitHub icon
        html.A(
            html.Img(
                src="https://cdn.simpleicons.org/github/ff8aff",
                style={
                    'height': '1.2vw',
                    'width': '1.2vw',
                    'marginRight': '0.5vw',
                    'verticalAlign': 'middle'
                }
            ),
            href=member["github"],
            target="_blank"
        ),

        # Email icon
        html.A(
            html.Img(
                src="https://cdn.simpleicons.org/gmail/ff8aff",
                style={
                    'height': '1.2vw',
                    'width': '1.2vw',
                    'verticalAlign': 'middle'
                }
            ),
            href=member["email"],
            target="_blank"
        )
    ], style={
        'marginBottom': '0.8vw',
        'display': 'flex',
        'alignItems': 'center'
    })

layout = html.Div(
    style={
        'minHeight': '100vh',
        'background': 'radial-gradient(ellipse at 70% 70%, #b100b1 0%, #1a001a 80%)',
        'color': '#fff',
        'padding': '0',
        'margin': '0',
        'display': 'flex',
        'alignItems': 'center',
        'justifyContent': 'flex-start',
        'fontFamily': "'DM Mono', monospace"
    },
    children=[
        html.Div(
            style={
                'marginLeft': '7vw',
                'maxWidth': '800px',
                'display': 'flex',
                'flexDirection': 'column',
                'justifyContent': 'center',
                'height': '100vh'
            },
            children=[
                # TITLE
                html.H1(
                    "SINGLE CELL VISUALIZATION SUITE",
                    style={
                        'fontFamily': "'Montserrat', sans-serif",
                        'fontWeight': 200,
                        'fontSize': '5vw',
                        'letterSpacing': '0.14em',
                        'marginBottom': '2vw',
                        'lineHeight': '1.05',
                        'textTransform': 'uppercase',
                        'color': '#fff'
                    }
                ),
                # PROJECT INFO
                html.Div([
                    html.Div(
                        "Democratizing analytics with integrated visualization",
                        style={
                            'fontFamily': "'DM Mono', monospace",
                            'fontSize': '1.4vw',
                            'letterSpacing': '0.08em',
                            'marginBottom': '1vw',
                            'color': '#ccc',
                            'fontWeight': 400
                        }
                    ),
                    html.Div(
                        "An interactive dashboard for exploring and visualizing single-cell RNA-seq data. Includes dimensionality reduction, gene expression analysis, and more.",
                        style={
                            'fontFamily': "'DM Mono', monospace",
                            'fontSize': '1.2vw',
                            'letterSpacing': '0.03em',
                            'marginBottom': '2vw',
                            'color': '#ccc',
                            'fontWeight': 400
                        }
                    ),
                ]),

                # TEAM MEMBER(S) WITH ICONS
                html.Div(
                    [create_team_member_entry(member) for member in TEAM_INFO],
                    style={'marginBottom': '2vw', 'marginTop': '1vw'}
                ),

                # BUTTON
                dcc.Link(
                    html.Button(
                        'START EXPLORING   →',
                        style={
                            'background': 'transparent',
                            'color': '#ff8aff',
                            'border': '2px solid #ff8aff',
                            'padding': '1.2vw 2.5vw',
                            'fontSize': '1.5vw',
                            'borderRadius': '6px',
                            'marginTop': '2vw',
                            'cursor': 'pointer',
                            'letterSpacing': '0.12em',
                            'fontWeight': 400,
                            'fontFamily': "'DM Mono', monospace",
                            'transition': 'background 0.2s, color 0.2s'
                        }
                    ),
                    href='/dashboard'
                )
            ]
        )
    ]
)