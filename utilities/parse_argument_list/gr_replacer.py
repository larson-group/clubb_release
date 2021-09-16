#!/usr/bin/env python3
# @author: Keegan Rhodes
'''
Biggest issue with this currently is
the fact that private methods are being overwritten if they are the same name
Challenges:


    functions with only one parameter having gr being inserted with gr,
    dupilcate function names

    Certain subourtines/functions being marked as changed when they weren't and
    other being changed without being marked by appearing in the list of changed

    If wants to be used more widescale more testing is defientley going to be needed


    Add feature that checks one line above to see how much preceeding whitespace
'''

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


def being_declared(string):
    """
    Helper method used to see if the function or subroutine is being defined
    :param string: the string being checked against the forbidden words
    :return: a boolean indicating if it is being declared or called
    """

    if 'write' in string or 'function' in string or 'subroutine' in string \
            or 'character' in string or 'if' in string or 'result' in string:
        being_declared = True
    else:
        being_declared = False
    return being_declared


def already_visited(string):
    """
    Helper method used to identify if a subroutine call or definition has
    already been visited by the script in another instance
    :param string: The call or definition of a subroutine/function
    :return: a boolean indicating if it has been visited already or not
    """
    separated = string.partition('(')[2]
    if separated.replace(' ', '').replace('(', '')[:2] == 'gr':
        visited = True
    else:
        visited = False

    return visited


def clean_list(list_of_subroutines):
    """
    Simple function that just removes empty strings from the lsit of subroutines
    :param list_of_subroutines: the list of subroutines and functions affected
                                by adding gr into the decleration
    :return: a cleaner version of the list of subroutines without empty strings
    """
    for word in list_of_subroutines:
        if word == '':
            list_of_subroutines.remove(word)
    return list_of_subroutines


def get_subroutine_name(string):
    """
    Function that splits a string on everything before the first (
    This may cause issues if there are more than one ( however, this is
    intended to be used on a subroutine/function decleration statement
    so there should only be one ( character in the entire string.
    :param string: the string which contains the name of a subroutine/function
                   which is wanting to be extracted.
    :param name: if this method is being called on a function or a
                 subroutine.
    :return: returns the name of the subroutine with no extra characters
    """
    separated = string.partition('(')
    separated = separated[0].replace('subroutine', ''). \
        replace(',', '').replace('&', '').replace('call', ''). \
        replace('pure', '').replace('elemental', '').replace('function', ''). \
        replace(' ', '').replace('\n', '')

    return separated


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


def made_to_grid_class(content, index):
    """
    Checks to see if the index provided has made it to grid_class or not
    :param content: the content containing all aspects of the file being read
    :param index: the index indicating which line to look at
    :return: a boolean indicating if grid class has been reached or not
    """
    grid_class_index = content[index].find('clubb_api_module')
    grid_class_holder = content[index][grid_class_index + 10:]
    gotten_to_use_grid_class = False
    if 'gr' in grid_class_holder and not is_comment(content[index], 'gr'):
        gotten_to_use_grid_class = True
    elif 'gr' not in grid_class_holder:
        j = index + 1
        while 'use' not in content[j]:
            if 'gr' in content[j] and not is_comment(content[j], 'gr'):

                testing_string = content[j].replace(' ', '').replace('\n', ''). \
                    partition('!')[0]
                if testing_string == 'gr,&' or testing_string == 'gr' or \
                        testing_string == 'grid' or testing_string == 'grid,&':
                    gotten_to_use_grid_class = True
                    break
            else:
                gotten_to_use_grid_class = False
            j += 1
    return gotten_to_use_grid_class


def add_gr_to_arg_list(content, list_of_subroutines):
    """
    This method will locate a subroutine that needs to have gr put into its
    argument list and then add it in. Additionally, if a change is made to
    a subroutine's arg list, its name will be flagged and put into
    list_of_subroutines for later modification of calls.
    :param content: the content of the file being edited
    :param list_of_subroutines: list of subroutines/functions being modified
    """
    name_index = -1
    content_index = -1
    gotten_to_use_grid_class = False
    new_subroutine_found = False
    for i in range(len(content)):
        if "subroutine" in content[i] and '(' in content[i] \
                and not is_comment(content[i], 'subroutine'):
            name_index = i
            new_subroutine_found = True
            content_index = i
        elif 'subroutine' in content[i] and '(' in content[i + 1] \
                and not is_comment(content[i], 'subroutine'):
            name_index = i
            new_subroutine_found = True
            content_index = i + 1
        elif "function" in content[i] and '(' in content[i] \
                and not is_comment(content[i], 'function'):
            name_index = i
            new_subroutine_found = True
            content_index = i
        elif 'function' in content[i] and '(' in content[i + 1] \
                and not is_comment(content[i], 'function'):
            name_index = i
            new_subroutine_found = True
            content_index = i + 1

        if 'end subroutine' in content[i] or 'end function' in content[i]:
            new_subroutine_found = False

        if 'use grid_class' in content[i] and new_subroutine_found:
            gotten_to_use_grid_class = made_to_grid_class(content, i)

        if gotten_to_use_grid_class and not already_visited(content[content_index]):

            name = get_subroutine_name(content[name_index])
            insert_gr_to_list(content, content_index, '(')

            if name not in list_of_subroutines:
                list_of_subroutines.append(name)
            gotten_to_use_grid_class = False
            new_subroutine_found = False
        elif already_visited(content[content_index]):
            gotten_to_use_grid_class = False
            new_subroutine_found = False


def add_grid_typing(content):
    """
    This method will add the grid typing string into the function/subroutine
    so that way everything will be in scope for when gr is passed in
    :param content: the content of the file being edited
    """
    gotten_to_implicit = False
    gotten_to_use_grid_class = False
    for i in range(len(content) - 1):
        if gotten_to_implicit and content[i + 1] != '    type (grid), target, intent(in) :: gr\n':
            content.insert(i, '\n    type (grid), target, intent(in) :: gr\n')
            gotten_to_implicit = False
            gotten_to_use_grid_class = False
        elif content[i + 1] == '    type (grid), target, intent(in) :: gr\n':
            gotten_to_implicit = False
            gotten_to_use_grid_class = False
        if 'use clubb_api_module' in content[i]:
            gotten_to_use_grid_class = made_to_grid_class(content, i)

        if gotten_to_use_grid_class and 'implicit none' in content[i]:
            gotten_to_implicit = True


def remove_gr(content):
    """
    This is the method that will remove 'use gr' from subroutines and replace
    it with 'use grid' instead
    :param content: the actual content of the file being edited
    """
    gotten_to_use_grid = False
    for i in range(len(content)):

        if gotten_to_use_grid:

            j = i
            while 'use' not in content[j] and not is_comment(content[j], 'gr') \
                    and 'implicit none' not in content[j]:
                testing_string = content[j].replace(' ', '').replace('\n', ''). \
                    partition('!')[0]
                if testing_string == 'gr,&':
                    content[j] = content[j].replace(content[j], '        grid, & ! Type\n')
                elif testing_string == 'gr':
                    content[j] = content[j].replace(content[j], '        grid ! Type\n')
                j += 1
            gotten_to_use_grid = False
        if 'use grid_class' in content[i]:
            if 'use grid_class, only: gr' in content[i] and not is_comment(content[i], 'use'):
                content[i] = content[i].replace(content[i], '    use grid_class, only: grid\n')
            gotten_to_use_grid = True


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


def prep_for_gr(file, path, list_of_subroutines):
    """
    Method which will call all other methods used to parse information
    :param file: the file being read/written to
    :param path: the path of the file
    :param list_of_subroutines: the running list of all subroutines being affected
    """
    clear_to_edit = False
    content = []
    file_to_edit = open(path + '/' + file, 'r')
    lines = file_to_edit.readlines()
    for word in lines:
        content.append(word)
    file_to_edit.close()
    for word in content:  # This whole block needs to be removed somehow without breaking
        if "use grid_class" in word:
            clear_to_edit = True
            break
    if clear_to_edit:
        add_gr_to_arg_list(content, list_of_subroutines)
        clean_list(list_of_subroutines)
        remove_gr(content)
        add_grid_typing(content)
    save_content(path, file, content)
    content.clear()


def insert_gr_to_list(content, index, name):
    """
    Helper method used to add gr to the calling list
    :param content: the content of a given line
    :param index: the index for which line to look at in content
    :param name: the name of the function/subroutine
    """
    name_index = content[index].find(name)

    while content[index][name_index + len(name)] == ' ':
        name_index += 1
    if content[index][name_index + len(name)] == '(':
        name_index += 1

    first_part = content[index][:(name_index + len(name))]

    second_part = content[index][(name_index + len(name)):]
    counter = 1
    while second_part[0] == ' ':
        
        second_part = content[index][(name_index + len(name) + counter): ]
        counter += 1
    if first_part[len(first_part)-1] == ' ':
        first_part += 'gr, '
    else:
        first_part += ' gr, '
    first_part += second_part
    content[index] = first_part


def add_gr_to_call(file, path, list_of_subroutines):
    """
    Method used to add gr to every calling of a subroutine or function
    :param file: the file being read/written to
    :param path: the path of the file
    :param list_of_subroutines: running list of subroutines affected by change
    """
    content = []

    file_to_edit = open(path + '/' + file, 'r')
    lines = file_to_edit.readlines()
    for word in lines:
        content.append(word)
    file_to_edit.close()
    for i in range(len(content) - 1):  # TODO definitley refactor this is horrible
        if first_non_whitespace_char(content[i]) != '!':
            test_string = content[i]
            if '=' in test_string:
                test_string = test_string.partition('=')[2]
                if '+' in test_string:
                    test_string = test_string.partition('+')[2]

            if 'call' in content[i] and not already_visited(content[i]) \
                    and not already_visited(content[i + 1]):  # TODO refactor

                name = get_subroutine_name(content[i])
                if name in list_of_subroutines and '(' in content[i]:
                    insert_gr_to_list(content, i, name)
                elif name in list_of_subroutines and '(' in content[i + 1]:
                    insert_gr_to_list(content, i + 1, '(')
            elif not already_visited(test_string) and not already_visited(content[i + 1]):  # TODO refactor

                for function_name in list_of_subroutines:
                    if function_name in content[i] and '(' in content[i] and not \
                            being_declared(content[i]) and not is_comment(content[i], '('):
                        insert_gr_to_list(content, i, function_name)

    save_content(path, file, content)


def add_gr_retroactively(file, path, list_of_subroutines):
    changed = False
    unedited_content = []
    content = []
    file_to_edit = open(path + '/' + file, 'r')
    lines = file_to_edit.readlines()
    for line in lines:
        content.append(line)
        unedited_content.append(line)
    subroutine_located = False
    content_index = -1
    for i in range(len(content)):
        if 'subroutine' in content[i] and not is_comment(content[i], 'subroutine'):
            subroutine_located = True
            name = get_subroutine_name(content[i])

            content_index = i
        elif 'function' in content[i] and not is_comment(content[i], 'function'):
            subroutine_located = True
            name = get_subroutine_name(content[i])
            content_index = i
        if 'end subroutine' in content[i] or 'end function' in content[i]:
            subroutine_located = False
            content_index = -1

        if subroutine_located and 'gr,' in content[i] \
                and not is_comment(content[i], 'gr,') and name not in list_of_subroutines:

            while 'use' not in content[content_index] and not is_comment(content[content_index], 'use'):
                content_index += 1
            while content[content_index] != '\n':
                content_index += 1
            content.insert(content_index, '    use grid_class, only: grid\n')
            subroutine_located = False

    if unedited_content != content:
        save_content(path, file, content)
        changed = True
    return changed


def main():
    changes_made = False
    list_of_subroutines = []
    previous_file_if_applicable = input("If previously ran on another directory\n\
enther the path of the functions_affected.txt.file. Leave blank otherwise. ")
    path = '/home/rhodesk/Desktop/Projects/clubb/src'
    if previous_file_if_applicable != '':
        f = open(previous_file_if_applicable + '/functions_affected.txt', 'r')
        for line in f:
            list_of_subroutines.append(line.replace('\n', ''))

    files = [f for f in listdir(path) if isfile(join(path, f))]
    for file in files:
        prep_for_gr(file, path, list_of_subroutines)

    for file in files:
        add_gr_to_call(file, path, list_of_subroutines)

    for file in files:
        check = add_gr_retroactively(file, path, list_of_subroutines)
        if check:
            changes_made = True

    while changes_made:
        changes_made = False
        for file in files:
            prep_for_gr(file, path, list_of_subroutines)
        for file in files:
            add_gr_to_call(file, path, list_of_subroutines)
        for file in files:
            check = add_gr_retroactively(file, path, list_of_subroutines)
            if check:
                changes_made = True

    f = open(path + '/functions_affected.txt', 'a')
    # Line used to output all functions thought to be affected by changes
    for name in list_of_subroutines:
        f.write(name + '\n')
    f.close()


# FILES THAT BREAK IT
# numerical_check.F90
# advance_xm has two use grid class statements for some reason
# advance_windm having multiple gr's being thrown in
# gr being thrown into the middle of a variable name


if __name__ == "__main__":
    main()

