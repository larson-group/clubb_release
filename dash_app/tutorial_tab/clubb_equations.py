"""CLUBB-equations tutorial page."""

from dash import html

from .clubb_equations_demo.app import build_layout as build_equation_layout


def build_page():
    """Embed the interactive equation quick reference in the Tutorial tab."""
    return html.Div(build_equation_layout(), className="tutorial-equations-page")
