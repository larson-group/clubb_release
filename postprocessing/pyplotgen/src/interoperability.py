# -*- coding: utf-8 -*-
"""
:Author: Steffen Domke
:Date: March 2020

This file contains functions used to guarantee operability on
different operating systems.
"""
import os
import re

def clean_path(string):
    """
    Removes characters from a path string
    that are not valid for the operating system

    :param string: Path string to have characters removed
    :return: a cleaned version of the path
    """
    string = string.replace('.', '')
    string = string.replace(',', '')
    string = string.replace(' ','_')
    string = string.replace('*','x')
    if 'win' in os.name.lower() or 'nt' in os.name.lower():
        # Replace all ':' except the one specifying the drive
        # (Uses look ahead and look back patterns)
        string = re.sub(r'(?<![A-Z]):(?!\\)', '-', string)
    else:
        string = string.replace(':', '-')
    return string