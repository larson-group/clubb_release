# -*- coding: utf-8 -*-
"""
:Author: Steffen Domke
:Date: March 2020

This file contains functions used to guarantee operability on
different operating systems.
"""
import os
import re

def clean_path(input_string):
    """
    Removes characters from a path string
    that are not valid for the operating system

    :param input_string: Path string to have characters removed
    :return: A cleaned version of the path
    """
    # string = string.replace('.', '')
    input_string = input_string.replace(',', '')
    input_string = input_string.replace(' ', '_')
    input_string = input_string.replace('*', 'x')
    if 'win' in os.name.lower() or 'nt' in os.name.lower():
        # Replace all ':' except the one specifying the drive
        # (Uses look ahead and look back patterns)
        input_string = re.sub(r'(?<![A-Z]):(?!\\)', '-', input_string)
    else:
        input_string = input_string.replace(':', '-')
    return input_string

def clean_title(string):
    """
    Cleans plot titles so they can be used as part of a file name.
    Matplotlib plot titles may contain latex code enclosed in '$' characters.
    This code should be removed entirely to avoid having invalid characters in file names.
    
    :param string: Title of a plot which is used to generate a file name
    :return: A cleaned version of the input string
    """
    return re.sub('\$[^$]*\$','',string)