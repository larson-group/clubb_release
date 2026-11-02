"""Provide dashboard colors for Plotly figures.

Call figure_theme with the selected theme when building a figure. It returns
layout settings for the canvas, font and axes; trace colors remain the
responsibility of the figure so their scientific meanings are preserved.
"""

def figure_theme(theme):
    """Return Plotly layout colors for the selected dashboard theme."""
    if theme == 'light':
        return dict(template='plotly_white', paper_bgcolor='#ffffff', plot_bgcolor='#ffffff',
                    font=dict(family='Arial, sans-serif', color='#172d3e'),
                    xaxis=dict(gridcolor='#dce7ed', zerolinecolor='#c7d6e0'),
                    yaxis=dict(gridcolor='#dce7ed', zerolinecolor='#c7d6e0'))
    return dict(template='plotly_dark', paper_bgcolor='#151f30', plot_bgcolor='#151f30',
                font=dict(family='Arial, sans-serif', color='#edf3fa'),
                xaxis=dict(gridcolor='#35465e', zerolinecolor='#536980'),
                yaxis=dict(gridcolor='#35465e', zerolinecolor='#536980'))
