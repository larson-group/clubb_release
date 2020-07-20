# -*- coding: utf-8 -*-
'''
Child class derived from Panel.
Takes the same variable definitions as a Panel for profile plots.
But needs the unaveraged netcdf data.
The changed plot routine outputs animations of profile plots instead of simple profile plots.

:author: Steffen Domke
:date: July 2020
'''
import os
from datetime import datetime
import matplotlib.pyplot as plt
from Panel import Panel

from config import Style_definitions
from src.interoperability import clean_path, clean_title

class AnimationPanel(Panel):
    """
    ContourPanel class derived from Panel
    The difference is that the plot routine will generate a timeheight plot instead of a standard profile plot.
    """

    def __init__(bla):
        """
        Init, same as Panel

        :param bla: 
        :return: 
        """
        super(bla)

    def plot(bla):
        """
        New version of plot routine to generate movies of profiles.

        :param bla: 
        :return: None
        """
        pass