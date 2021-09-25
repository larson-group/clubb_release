#!/usr/bin/env python3

from os import listdir
from os.path import isfile, join




def is_comment(string, delimiter):
    """
    Helper method to determine if something is a comment or not
    :param string: the string being checked to see if it's a comment
    :param delimiter: a string for what you what to check to see if it comes
                    before or after the '!'
    :return: a boolean indicating if it is a comment
    """
    testing_index = string.find('!')
    green_light_index = string.find(delimiter)
    if testing_index == -1 or green_light_index < testing_index:
        comment = False
    else:
        comment = True
    return comment



def first_non_whitespace_char(string):
    """
    Simple function that checks to see what the first non whitespace char is
    :param string: the string that is going to be checked
    :return: returns the first non whitespace char of string
    """
    char = None
    for i in range(len(string)):
        if string[i] != ' ':
            char = string[i]
            break
    return char

def push_all_check(string):
    check = False
    forbidden = ('+', '-', '*', '/', '=', 'for', 'do', 'if', 'call')
    for word in forbidden:
        if word in string and not is_comment(string, word):
            if word == 'for' or word == 'do' or word == 'if' or word == 'call':
                if use_check(string, word):
                    check = True
            else:    
                check = True
    return check

def enforce_indent(file, path):
    push_all = False
    content = []
    gotten_to_subroutine = False
    file_to_edit = open(path + '/' + file, 'r')
    lines = file_to_edit.readlines()
    for line in lines:
        content.append(line)
    for i in range(len(content)):
        if push_all_check(content[i]) and '::' not in content[i]:
            push_all = False
        if gotten_to_subroutine:
            num_whitespace = len(content[i]) - len(content[i].strip(' '))
            if ('use' in content[i] or 'implicit' in content[i]) and first_non_whitespace_char(content[i]) != '!'\
            and num_whitespace < 4:
                push_all = True
            if push_all and 'return' not in content[i] and 'end' not in content[i]\
            and '#' not in content[i] and content[i].replace(' ', '').replace('\n', '')\
            != '' and first_non_whitespace_char(content[i]) != '!':
                first_part = '  '
                first_part += content[i]
                content[i] = first_part
                
            if 'use' in content[i] and not is_comment(content[i], 'use') and\
            use_check(content[i], 'use'):
                j = i + 1
                while 'use' not in content[j] and 'implicit' not in content[j]:
                    space = len(content[j]) - len(content[j].strip(' '))
                    if '#' not in content[j] and content[j].replace(' ', '').\
                    replace('\n', '') != '' and space < 8 and\
                    first_non_whitespace_char(content[j]) != '!':
                        first_part = '  '
                        first_part += content[j]
                        content[j] = first_part
                    j += 1
                    

       
        if 'subroutine' in content[i] and not is_comment(content[i], 'subroutine'):
            gotten_to_subroutine = True
        elif 'function' in content[i] and not is_comment(content[i], 'function'):
            gotten_to_subroutine = True
        
        if 'end subroutine' in content[i] or 'end function' in content[i]:
            gotten_to_subroutine = False
            push_all = False
        
    save_content(path, file, content)      

def use_check(string, delimiter):
    real = False
    string = string.replace(' ', '')
    index = string.find(delimiter)
    if index == 0:
        real = True
    return real
    
        
def save_content(path, file, content):
    """
    How changes are actually made to the underlying file without
    editing the file in place
    :param path: the absolute path of the file
    :param file: the name of the file
    :param content: the structure which houses all the lines of the file
    """
    file_to_edit = open(path + '/' + file, 'w')
    for word in content:
        file_to_edit.write(word)
    file_to_edit.close()
    
    
def main():
    
    path = '/home/rhodesk/Desktop/Projects/changed/clubb/src/CLUBB_core'
    files = [f for f in listdir(path) if isfile(join(path, f))]
    for file in files:
        enforce_indent(file, path)
    



if __name__ == '__main__':
    main()