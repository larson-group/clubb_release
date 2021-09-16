#!/usr/bin/env python3
# @author: Keegan Rhodes

'''
Plan of attack right now is to make a script which will identify all subroutines
whose arg list will need to be updated either directly or indirectly and then
generate a list based on that information. Later, that list will be used to
update via another script.
'''

from os import listdir
from os.path import isfile, join


def remove_comments(string):
    comment_index = string.find('!')
    if comment_index != -1:    
        string = string[:comment_index]
    return string

def remove_whitespace(string):
    string = string.lstrip()
    string = string.rstrip()
    string = ' '.join(string.split())
    return string

def purge(string):
    if 'use' not in string and 'subroutine' not in string and 'call' not in string\
    and ', &' not in string and 'implicit' not in string:
        string = ''
    return string

def prep_data(content):
    
    content = list(map(remove_comments, content))
    
    content = list(map(remove_whitespace, content))
    content = list(map(purge, content))
    content = list(map(lambda string: string.replace('&', '').replace(',', ''), content))
    content = list(map(remove_whitespace, content))

    content = list(filter(None, content))
    return content


def spot_check(string):
    stats = ''
    checked = False
    targets = ('stats_zm', 'stats_zt', 'stats_sfc', 'stats_lh_zt', 'stats_lh_sfc', 'stats_rad_zt', 'stats_rad_zm')
    for i in range(len(targets)):
        if string == targets[i]:
            checked = True
            stats = targets[i]
    return checked, stats


def generate_list(list_affected, content):
    stat_names = []
    name = ''
    for i in range(len(content)):
        if 'subroutine' in content[i]:
            name = content[i].replace('subroutine', '').replace(' ', '')\
            .replace('pure', '').replace('elemental', '')
            if name.find('(') != -1:
                index = name.find('(')
            elif name.find('&') != -1:
                index = name.find('&')
            else:
                index = len(name)
            name = name[:index]
        checked, stats = spot_check(content[i])
        if stats != '' and stats not in stat_names:
            stat_names.append(stats)
        if len(stat_names) != 0:
            if 'end' in name and name not in list_affected:
                name = name.replace('end', '')
                list_affected[name] = stat_names
             
                stat_names = []


def retroactive_case(list_affected, file, path):
    content = []
    file_to_edit = open(path + '/' + file, 'r')
    lines = file_to_edit.readlines()
    for line in lines:
        content.append(line)
    file_to_edit.close()
    modded = False
    subroutine_name = ''
    name = ''
    content = prep_data(content)
    for i in range(len(content)):
        if 'subroutine' in content[i] and 'end' not in content[i]:
            subroutine_name = content[i].replace('subroutine', '').\
            replace('pure', '').replace('elemental', '').replace(' ', '')
            if subroutine_name.find('(') != -1:
                index = subroutine_name.find('(')
            elif subroutine_name.find('&') != -1:
                index = subroutine_name.find('&')
            else:
                index = len(subroutine_name)
            subroutine_name = subroutine_name[:index]
        if 'call' in content[i]:
            name = content[i].replace('call', '')
            name = remove_whitespace(name)
            if name.find(' ') != -1:
                index = name.find(' ')
            elif name.find('(') != -1:
                index = name.find('(')
            elif name.find('&') != -1:
                index = name.find('&')
            else:
                index = len(name)
            name = name[:index]
            if name in list_affected and subroutine_name not in list_affected:
                list_affected[subroutine_name] = list_affected[name]
                modded = True
    return modded
            



def locate_subroutines(list_affected, file, path):
    content = []
    file_to_edit = open(path + '/' + file, 'r')
    lines = file_to_edit.readlines()
    for line in lines:
        content.append(line)
    file_to_edit.close()
    
    content = prep_data(content)
    generate_list(list_affected, content)
    


                                                              
def main():
    list_affected = {}
    path = '/home/rhodesk/Desktop/Scripts/unsure_output'
    files = [f for f in listdir(path) if isfile(join(path, f))]
    for file in files:
        locate_subroutines(list_affected, file, path)
    modded = True
    while modded:
        for file in files:
            modded = retroactive_case(list_affected, file, path)
 
    f = open(path + '/the_affected.txt', 'w')
    for key in list_affected:
        f.write(key + ' ' + str(list_affected[key]) + '\n')
    f.close()
    
    


if __name__ == '__main__':
    main()