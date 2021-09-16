#!/usr/bin/env python3
# @author: Keegan Rhodes


from os import listdir
from os.path import isfile, join
from re import sub
import copy




'''
Known issues:
    Issue with comment throwing off the reverse tehcnique
    Issue with subroutines that don't have () in their calls or definitions
'''

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
    pattern = r'[0-9]'
    targets = ('stats_zm', 'stats_zt', 'stats_sfc', 'stats_lh_zt', 'stats_lh_sfc', 'stats_rad_zt', 'stats_rad_zm')
    if 'use' not in string and 'subroutine' not in string and 'call' not in string\
    and ', &' not in string and 'implicit' not in string and sub(pattern, '', string).strip() not in targets:
        string = ''
    return string

def prep_data(content):
    
    content = list(map(remove_comments, content))
    content = list(map(remove_whitespace, content))
    content = list(map(purge, content))
    # content = list(map(lambda string: string.replace('&', '').replace(',', ''), content))
    content = list(map(remove_whitespace, content))

    content = list(filter(None, content))
    return content



def add_the_stuff(list_affected, file, path):
    args_created = False
    name = ''
    content = []
    unchanged_content = []
    modded_content = {}
    file_to_edit = open(path + '/' + file, 'r')
    lines = file_to_edit.readlines()
    count = 0
    
    #Adds line numbers and junk
    for line in lines:
        unchanged_content.append(line)
        string = str(count) + ' ' + line
        content.append(string)
        count += 1
    file_to_edit.close()
    
    content = prep_data(content)
   
    #Block that flags which line needs to be modified 
    for i in range(len(content)):
        if ('subroutine' in content[i] and 'end subroutine' not in content[i])\
            or 'call' in content[i]:
            name = content[i].replace(' ', '').replace('call', '')\
            .replace('subroutine', '').replace('pure', '').replace('elemental', '')

            index = 0
            while not name[index].isalpha():
                index += 1
            name = name[index:]
            
            
            if name.find('(') != -1:
                index = name.find('(')
            elif name.find('&') != -1:
                index = name.find('&')
            else:
                index = len(name)
            name = name[:index]

            if name in list_affected:
               
                index = content[i].find(' ')
                line_number = content[i][:index]
                
                modded_content[line_number] = list_affected[name]
                
        
    #Block that will add the needed variables in arg lists and calls to
    for i in range(len(unchanged_content)):
        if str(i) in modded_content:
            content_index = i

            if parameter_checker(unchanged_content[content_index]) and '&' not in unchanged_content[content_index]:
                temp_string = unchanged_content[content_index][:len(unchanged_content[content_index])-1]
                if '( )' not in temp_string:
                    temp_string += ' ( )\n'
                unchanged_content[content_index] = temp_string
                args_created = True
            
            elif '(' not in unchanged_content[content_index]:
                content_index += 1

            temp_list = modded_content[str(i)]
          
            
            
            num_spaces = len(unchanged_content[content_index+1]) - len(unchanged_content[content_index+1].lstrip())
            
            
            
            temp_string = unchanged_content[content_index]
            
            if ')' in unchanged_content[content_index] and not args_created\
            and not is_comment(unchanged_content[content_index], ')')\
            and not is_array(unchanged_content[content_index], ')'):
                
                index = remove_comments(temp_string).rfind(')')-1
                temp_string = temp_string[:index]
                temp_string += ', & !intent(in)\n'
            
            
            elif ')' in unchanged_content[content_index] and args_created\
            and not is_comment(unchanged_content[content_index], ')')\
            and not is_array(unchanged_content[content_index], ')'):
                
                index = remove_comments(temp_string).rfind(')')-1
                temp_string = temp_string[:index]
                num_spaces = 0
            
            for i in range(num_spaces-2):
                temp_string += ' '
                
            temp_string += ' '
            
            
            for j in range(len(temp_list)):
                
                if args_created:
                    temp_string += temp_list[j] + ','
                else:
                    temp_string += ' ' + temp_list[j] + ','
                
                args_created = False


            if ')' in unchanged_content[content_index]\
            and not is_comment(unchanged_content[content_index], ')')\
            and not is_array(unchanged_content[content_index], ')'):
                
                index = temp_string.rfind(',')
                temp_string = temp_string[:index]
                temp_string += ' ) ! intent(inout)'

            else:
                temp_string += ' & ! intent(inout)'
            
            unchanged_content[content_index] = temp_string + '\n'
    save_content(path, file, unchanged_content)





def remove_the_things(list_affected, file, path):
    file_to_edit = open(path + '/' + file, 'r')
    lines = file_to_edit.readlines()
    modded_content = {}
    content = []
    unchanged_content = []
    count = 0
    for line in lines:
        string = str(count) + ' ' + line
        content.append(string)
        unchanged_content.append(line)
        count += 1
    file_to_edit.close()
    
    content = prep_data(content)
    

    for i in range(len(content)):
        checked, stats = spot_check(content[i])
        if checked:
            temp = copy.deepcopy(content[i][:content[i].find(' ')])
            modded_content[temp] = ''
            # if ',' not in content[i] and 'use' not in content[i-1]:
            #     while content[i-temp_count].replace(' ', '').replace('\n', '') != '':
            #         temp_count += 1
            #     
            #     temp_string = content[i-temp_count][content[i-temp_count].find(' ')+1:]
            #     temp_string = temp_string[:temp_string.rfind(',')]
            #     temp = copy.deepcopy(content[i][:content[i].find(' ')])
            #     modded_content[temp] = '        ' + temp_string + '\n'
            #     content[i-1] = '        ' + temp_string +'\n'
            # elif ',' not in content[i] and 'use' in content[i-1]:
            #     temp_string = content[i-1][content[i-1].find(' ')+1:]
            #     temp_string = temp_string[:temp_string.rfind(',')]
            #     temp = copy.deepcopy(content[i][:content[i].find(' ')])
            #     modded_content[temp] = temp_string
            #     content[i-1] = ''
            # content[i] = ''
            
    for key in modded_content:
        unchanged_content[int(key)] = modded_content[key]
    save_content(path, file, unchanged_content)
        
          
            
def add_the_other_stuff(list_affected, file, path):
    subroutine_name = ''
    file_to_edit = open(path + '/' + file, 'r')
    passed_subroutine = False
    lines = file_to_edit.readlines()
    content = []
    unchanged_content = []
    modded_content = {}
    count = 0
    for line in lines:
        string = str(count) + ' ' + line
        content.append(string)
        unchanged_content.append(line)
        count += 1
        
    content = prep_data(content)
 
    for i in range(len(content)):
        if 'subroutine' in content[i] and 'end subroutine' not in content[i]:
            subroutine_name = content[i][content[i].find(' ')+1:].replace(' ', '').\
            replace('subroutine', '').replace('pure', '').replace('elemental', '')
            if subroutine_name.find('(') != -1:
                index = subroutine_name.find('(')
            elif subroutine_name.find('&') != -1:
                index = subroutine_name.find('&')
            else:
                index = len(subroutine_name)
            subroutine_name = subroutine_name[:index]
            if subroutine_name in list_affected:
                passed_subroutine = True
        
        if 'implicit none' in content[i] and passed_subroutine:

            index = content[i][:content[i].find(' ')]

            modded_content[int(index)-1] = '\n    use stats_type, only: stats ! Type\n\n'
            temp = '\n    type (stats), target, intent(inout) :: &\n'
            for j in range(len(list_affected[subroutine_name])):
                if j != len(list_affected[subroutine_name])-1:
                    temp += '      ' + list_affected[subroutine_name][j] + ', &\n'
                elif j == len(list_affected[subroutine_name])-1:
                    temp += '      ' + list_affected[subroutine_name][j] + '\n\n'
            modded_content[int(index)+1] = temp
            passed_subroutine = False
            
    
    for key in modded_content:
        unchanged_content[int(key)] = modded_content[key]
        
    save_content(path, file, unchanged_content)
    
    

def spot_check(string):
    string = string.replace(',' , '').replace('&', '').strip()
    stats = ''
    string = string[string.find(' ')+1:]
    checked = False
    targets = ('stats_zm', 'stats_zt', 'stats_sfc', 'stats_lh_zt', 'stats_lh_sfc', 'stats_rad_zt', 'stats_rad_zm')
    for i in range(len(targets)):
        if string == targets[i]:
            checked = True
            stats = targets[i]
    return checked, stats




def is_array(string, delimiter):
    string = string.replace(' ', '')
    index = string.find(delimiter)
    is_array = False
    try:
        if string[index+1] != '!' and string[index+1] != '\n':
            is_array = True
    except:
        print('there was an error that had occured when parsing')
    return is_array
    



def is_comment(string, delimiter, method='f'):
    testing_index = string.find('!')
    if method == 'r':
        green_light_index = string.rfind(delimiter)
    else:    
        green_light_index = string.find(delimiter)
    if testing_index == -1 or green_light_index < testing_index:
        comment = False
    else:
        comment = True
    return comment



def parameter_checker(string):
    create_args = False
    temp_string = ''
    index = string.find('(')
    if index == -1:
        return True
    while string[index] != '\n':
        temp_string += string[index]
        index += 1

    if temp_string.replace(' ', '').replace('\n', '') == '()':
        create_args = True
    return create_args
                    
                                    
def save_content(path, file, content):
    file_to_edit = open(path + '/' + file, 'w')
    for word in content:
        file_to_edit.write(word)
    file_to_edit.close()



def main():
    list_affected = {}
    try:
        affected_path = '/home/rhodesk/Desktop'
        affected_file = open(affected_path + '/the_affected.txt', 'r')
        lines = affected_file.readlines()
        for line in lines:
            space_index = line.find(' ')
            name = line[:space_index]
            temp_string = line[space_index+1:].replace('\'', '').replace('[', '')\
            .replace(']', '').replace('\n', '').replace(',', '')
            temp_list = temp_string.split(' ')
            list_affected[name] = temp_list
            
    except FileNotFoundError:
        print('Could not locate the file specified, exiting')
        
    path = '/home/rhodesk/Desktop/Projects/clubb/src'
    files = [f for f in listdir(path) if isfile(join(path, f))]
    for file in files:
        remove_the_things(list_affected, file, path)
    for file in files:
        add_the_stuff(list_affected, file, path)
    for file in files:
        add_the_other_stuff(list_affected, file, path)


if __name__ == '__main__':
    main()