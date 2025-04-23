#!/usr/bin/env python3

# Script replacing all integer err_code variables with err_info_type struct
# Replacement changes have been committed to master
# TODO: (Parts that are missing/not working)
# 1. Adding the err_header to all error messages
# 2. I had to re-do the initialization scheme since I had no complete overview about the handling at top level (standalone, clubb tuner, etc.)
# 3. Distinction between column-internal and global err_header

import sys
import os
import re
import argparse

clubb_src_path = "/home/sdomke/workspace/clubb2/clubb/src"
use_pattern = ""
fortranFileExtensions = ('.f', '.f90')

#######################################################################################################################################
def main(path_list):
    path_list = sorted(path_list)
    print('path_list = ', path_list)
    f_list = []
    # path_list is expected to be a string or a list of strings.
    # If it's just a string wrap it into a list
    if not isinstance(path_list, list):
        path_list = [path_list]
    # Iterate through elements
    for path in path_list:
        # If path represents a Fortran source file, add it to the list
        if os.path.isfile(path) and path.lower().endswith(fortranFileExtensions):
            f_list.append(os.path.abspath(path))
        # If it's a directory. add all source files in directory tree with root path
        elif os.path.isdir(path):
            f_list.extend(find_all_src_files(path))

    print(f_list)
    print_all_file_endings(f_list)

    mod_files = []
    init_files = []
    for f in sorted(f_list):
        l_found, l_call_init = replace_all_err_code(f)
        if l_found:
            mod_files.append(f)
            print('File ', f, ' modified.')
        else:
            print('File ', f, ' NOT modified.')
        if l_call_init:
            init_files.append(f)

    print('Modified files:')
    print(sorted(mod_files))

    print('->Add init for err_info in these files:<-')
    print(sorted(init_files))

#end main

#######################################################################################################################################
def find_all_src_files(dir):
    """
    Returns a list of all Fortran source files in the file tree starting at dir

    This function is adapted from Stack Overflow:
    http://stackoverflow.com/questions/19932130/python-iterate-through-folders-then-subfolders-and-print-filenames-with-path-t

    :param dir: The directory to find the files in
    :return: A list of all the source file names
    """
    f_list = []
    # Walk through the directory tree
    for root, dirs, files in os.walk(dir):
        for fname in files:
            if fname.lower().endswith(fortranFileExtensions):
                f_list.append(os.path.join(root,fname))
    return f_list

#end find_all_src_files

#######################################################################################################################################
def print_all_file_endings(l):
    endings = set()
    for el in l:
        try:
            ending = el.rsplit('.', maxsplit=1)[1]
            endings.add(ending)
        except Exception as e:
            print('No . in ', el)
    print(endings)

#end print_all_file_endings

#######################################################################################################################################
def replace_all_err_code(fname):
    print('Replacing "err_code" in file ', fname)
    err_info_use = ['{}use err_info_type_module, only: &\n', '{}  err_info_type\t\t! Type\n', '\n']
    # Replace all mentions of err_code integer with err_info struct
    with open(fname, 'r') as f:
        content = f.readlines()
    # 1.a Replace err_code -> err_info in procedure signature
    # Build procedure pattern from blocks
    out = []
    # Signature flag: True  -> Trying to find err_code in signature
    #                 False -> Not in signature or err_code already replaced
    l_sig = False
    # Signature flag: True  -> Trying to find implicit none line to insert err_info use statement
    #                 False -> Not in procedure block or already replaced
    l_implicit_none = False
    l_mod1 = False
    l_mod2 = False
    l_mod3 = False
    l_mod4 = False
    l_mod5 = False

    for i, l in enumerate(content):
        if l_sig == False and l_implicit_none == False:
            out.append(l)
            if 'function' in l or 'subroutine' in l:
                # Procedure keyword found. Determine if it is a signature and if err_code is in parameter list
                line, end_i = combine_lines(i, content)
                if 'err_code' in line:
                    # Not in comment? Already excluded
                    print('"err_code" found in file '+ fname + ', lines ', i, '-', end_i, ': ' + line)
                    res = signature_check(line)
                    print(res)
                    if res:
                        print('Updating flags')
                        # Procedure with 
                        l_sig = True
                        l_implicit_none = True
                        # Written file will be modified
                        l_mod1 = True
        else:
            if l_sig == True:
                print('Trying to find "err_code" in this line:', l)
                # Still in signature, not yet found
                if re.match('[^!]*err_code',l):
                    print('"err_code" found -> Replacing')
                    # Word 'err_code' in signature found
                    # Replace 'err_code' with 'err_info'
                    print('old: ', l)
                    l = re.sub('err_code', 'err_info', l)
                    print('new: ', l)
                    # Append modified line
                    out.append(l)
                    # Mark as found in sig
                    l_sig = False
                else:
                    # Not yet found in signature
                    out.append(l)
            else:
                # 1.b If err_code in signature found:
                # Add err_info_type use: <Changed signature>...<Description comment>\n\n  use err_info_type, only : &\n    err_info_type
                print('Trying to find "implicit none" in this line:', l)
                # err_code in sig already replaced -> add use statement before implicit none
                res = re.search('^(\s*)implicit none', l)
                if res:
                    # Capture white space before 'implicit none'
                    prefix = res.groups()[0]
                    # Insert use statement, appropriately indented
                    out.extend([el.format(prefix) for el in err_info_use])
                    # Mark 'implicit none' as found
                    l_implicit_none = False

                # Add 'implicit none' line after either way
                out.append(l)

    if l_mod1:
        print('Err_code in signature was replaced')

    # 2. Replace variable definition:
    #    integer, intent(?) :: &...(err_code)... -> integer, intent(?) :: &...()...\n type(err_info_type)
    out, l_mod2, l_call_init = modify_var_def(out)

    if l_mod2:
        print('Err_code in variable definitions was replaced')

    # 3. Replace err_code in acc directives
    out, l_mod3 = replace_in_directives(out)

    if l_mod3:
        print('Err_code in directive was replaced')

    # 4. Replace any assignments: err_code = ? -> err_info%err_code = ?
    # 5. Replace any checks: if ( err_code == ? ) then ... -> if ( err_info%err_code == ? ) then ... 
    out, l_mod4 = replace_in_checks_and_assignments(out)

    if l_mod4:
        print('Err_code in checks or assignments was replaced')

    # 6. Add err_header to error messages -> Maybe only at top level?
    # 7. Replace in any calls: call(...err_code...) -> call(...err_info...)
    #                      or: out = func_name(...err_code...) -> func_name(...err_info...)
    out, lmod5 = modify_calls(out)

    if l_mod5:
        print('Err_code in procedure call was replaced')

    l_mod = l_mod1 or l_mod2 or l_mod3 or l_mod4 or l_mod5
    if l_mod:
        print('File was modified')
        with open(fname+'-mod', 'w') as outf:
            outf.writelines(out)
    print('Finished file ', fname)
    return l_mod, l_call_init

#end replace_all_err_code

#######################################################################################################################################
def combine_lines(n, lines):
    print('Calling combine lines on line ', n, ': ', lines[n])
    res = ''
    # "Skip" empty lines
    # 1st check: Line only contains white space and a comment
    # 2nd check: Line contains preprocessor directive
    # 3rd check: Line is empty (except for white space)
    if re.match(r"^\s*!(?!\$(acc|omp))", lines[n], re.IGNORECASE) or \
            lines[n].startswith("#") or \
            lines[n].strip() == "":
        print('Empty line replaced: ', lines[n])
        res = "&"                                        # Replace with '&' -> Combine with following lines
    else:
        res = lines[n]
    if '!$acc' in res or '!$omp' in res:
        res = res.strip()                                # Do not split in acc lines
    else:
        res = res.split("!")[0].strip()                  # Split off comments
    if res.endswith("&"):                                # If the line ends with a & -> continues in next line
        res = res.split("&")[0].rstrip()                 # Remove the &
        tmp, m = combine_lines(n + 1, lines)             # Call combine_lines on next line
        res += tmp                                       # Add the combined next line to this one
    else:
        m = n
    return res, m                                        # Return the reformatted line

#end combine_lines

#######################################################################################################################################
def signature_check(line):
    # Check given line if it contains a valid Fortran procedure definition
    # Define blocks
    prefixes = r'(\s*(pure|elemental))*'
    char_type = r'character\s*\(\s*len\s*=\s*\d+\s*\)'
    real_type = r'real\s*\(\s*kind\s*=\s*\w*\s*\)'
    int_type = r'integer'
    logical_type = r'logical'
    type_block = r'(' + char_type + r'|' + real_type + r'|' + int_type  + r'|' + logical_type + r')'
    dim_spec = r'\(\s*(:|\w+)\s*(,\s*(:|\w+)\s*)*\)'
    dim_key = r',\s*dimension\s*' + dim_spec
    name_pattern = r'[a-zA-Z]\w*'
    param_name_pattern = name_pattern + r'\s*(' + dim_spec + r')?'
    parameter_block = r'\(\s*' + param_name_pattern + r'(\s*,\s*' + param_name_pattern + r')*\s*\)'
    # Build procedure keyword blocks
    subroutine_block = r'subroutine'
    function_block = type_block + r'\s*' + dim_key + r'\s*function'
    # Combine blocks to procedure signature pattern up to and excluding the list of 
    signature_pattern = r'^' + prefixes + r'\s*(' + subroutine_block + r'|' + function_block + r')\s*' + name_pattern + r'\s*' + parameter_block

    return re.search(signature_pattern, line)

#end signature_check

#######################################################################################################################################
def modify_var_def(lines):
    # Find all variable definitions of err_code, remove, and replace with type(err_info_type)
    out = []
    # Marker indicating that we are inside a variable definition block containing a yet unreplaced 'err_code'
    l_find_err_code = False
    # Marker indicating that we are trying to find the right line to insert the new def block
    l_insert_def_block = False
    # Return value indicating if anything was changed in this procedure
    l_modified = False
    # Return value indicating if the init subroutine needs to be called in this file
    l_call_init = False
    prefix = ''
    postfix = ''
    # TODO: Append variable definition for err_info DONE
    # Define iterator so we can actually skip lines
    r = iter(range(len(lines)))
    for i in r:
        line = lines[i]
        # Default state: We're outside of a variable definition block
        if l_find_err_code == False and l_insert_def_block == False:
            if '::' in line:
                print('Var def block found', i)
                # Start of definition block found
                # Combine all lines and remove all comments
                combined_line, end_i = combine_lines(i, lines)
                print(i, end_i, combined_line)
                if 'err_code' in combined_line:
                    l_modified = True
                    print('"err_code" found', i)
                    # Isolate any other keywords for the definition so we can add it to the new block for err_info
                    keywords = re.search('^(.*?)integer(.*?)::', line)
                    try:
                        prefix = keywords.groups()[0]
                        postfix = keywords.groups()[1]
                        print(prefix)
                        print(postfix)
                        # TODO: check if intent(in) in prefix -> if not, init must be called (low prio, since does not happen often) DONE
                        if not re.search(r'intent\s*\(\s*in', postfix.lower()):
                            # err_code was not intent(out) -> notify to user to add init call in this file
                            print('init call for err_info needed after line ', i, ': ', line)
                            l_call_init = True
                    except Exception as e:
                        print('pre-/postfix pattern is not working', i)
                        sys.exit(1)
                    # Check if err_code is the only var defined in this block
                    only_var_check = re.search(r'::\s*err_code\s*$', combined_line)
                    if only_var_check or end_i == i:
                        print("Insert new block after this", i)
                        # If either err_code is the only variable in the block or the def block is only this line
                        # -> Append new def block directly after and we don't need to go looking for err_code after this
                        # err_code is only variable in this block -> Replace lines directly
                        if only_var_check is None:
                            # There are other vars -> Remove err_code, but keep line
                            print("Other vars -> Append modified line", i)
                            new_l = remove_var_in_line(line)
                            out.append(new_l)
                            # Add new empty line for nicer formatting
                            out.append('')
                        # else: Can just skip line(s)
                        # We want to skip the lines that are part of the continued line
                        # TODO: Check in output!!!
                        for n in range(i, end_i):
                            print("Skipping line", n)
                            next(r)
                        print("Inserting def block", i)
                        out.append(prefix + 'type(err_info_type)' + postfix + ':: &\n')
                        out.append(prefix + '  err_info\t\t! err_info struct containing err_code and err_header\n')
                    else:
                        print("Def block later", i)
                        # Other variables available and continued line
                        # -> New def block must be inserted later!
                        l_insert_def_block = True
                        # Check if same line or not
                        if 'err_code' in line:
                            print('err_code found in 1st line -> Append modified line', i)
                            # Check if 'err_code' is last var
                            # TODO: Preserve indentation?
                            # There is ALWAYS a comma after, due to checks
                            # Keep '&' in the same spot as before
                            start, end = re.search(r'err_code\s*,\s*', line).span()
                            new_l = re.sub(r'err_code\s*,\s*', ' '*(end-start), line)
                            # Add modified line either way
                            out.append(new_l)
                        else:
                            print('err_code in later line -> append old line and start search', i)
                            # Need to search for 'err_code' in next lines and remove
                            l_find_err_code = True
                            # Add old line
                            out.append(line)
                else:
                    out.append(line)
            else:
                # Not looking for err_code and not starting new def block -> Just copy/paste line
                out.append(line)
        else:
            # At least one of l_find_err_code or l_insert_def_block are True here
            # -> Either we want to find err_code in a def block or the end of the same def block to insert the new def block after
            # Isolate code part of line
            code_line = line.split('!', maxsplit=1)[0].strip()
            # TODO: We could have err_code in the last line -> both ifs and both checks are done
            # -> Be careful with variables and when to append!
            if l_find_err_code:
                # Search for line containing the definition of the var 'err_code'
                if 'err_code' in code_line:
                    print('err_code var def found', i)
                    print(line)
                    # Not searching for err_code anymore!
                    l_find_err_code = False
                    # Find indentation to add back in if necessary
                    indent = re.match(r'(\s*)[^\s]', line).groups()[0]
                    # Check if err_code is the only variable in this line
                    check = re.match(r'\s*err_code\s*(,\s*&)?$', code_line)
                    print(check)
                    if check:
                        print('err_code only var in this line -> skip', i)
                        # err_code is the only variable in this line
                        # If line does not end with '&', we must adjust the previous line
                        if check.groups() is None or check.groups()[0] is None or '&' not in check.groups()[0]:
                            print('Last line in def block -> Need to modify previous line', i)
                            print(out[-1])
                            # Last entry in variable list -> Skip this line and modify previous line
                            # Remove '&' from last line 
                            tmp = out[-1].replace('&', ' ')
                            print(tmp)
                            # Divide into code and comment
                            if '!' in tmp:
                                print('Splitting off comment part')
                                code, comment = tmp.rsplit('!', maxsplit=1)
                            else:
                                code = tmp
                                comment = ''
                            # Remove latest ',' in code
                            code = ' '.join(code.rsplit(',', maxsplit=1))
                            # Recombobulate code and comment parts
                            if comment:
                                out[-1] = code+'!'+comment
                            else:
                                out[-1] = code.rstrip()
                            print(out[-1])
                            # Add another line to separate definition blocks
                            out.append('\n')
                        else:
                            out.append('\n')
                            pass
                    else:
                        print('err_code is not the only var in this line -> Just remove', i)
                        new_l = remove_var_in_line(line)
                        out.append(new_l)
                else:
                    out.append(line)
            if l_insert_def_block:
                print('Trying to find end of integer block')
                print(code_line)
                # Trying to find end of integer def block
                if not code_line.endswith('&'):
                    print("Inserting def block", i)
                    out.append('\n')
                    out.append(prefix + 'type(err_info_type)' + postfix + ':: &\n')
                    out.append(prefix + '  err_info\t\t! err_info struct containing err_code and err_header\n')
                    l_insert_def_block = False

    return out, l_modified, l_call_init

#end modify_var_def

#######################################################################################################################################
def remove_var_in_line(target, var_name='err_code'):
    # Remove given var_name in target string, if it is NOY the only var in that line
    # TODO: Insert var_name into regular expressions
    split = re.match(r'(\s*)\w.*?(,)?\s*err_code\s*[,&]?', target)
    #                   ^ indentation
    #                        ^ first character of first var name
    #                             ^ Check if there is a ',' in front of err_code
    #                                                ^ Check if there is a ','or a '&' after err_code
    try:
        indent, pre_comma, post_char = split.groups()
        if pre_comma and post_char:
            print('Pre and post')
            # Limiting characters on both sides
            new_l = re.sub(r',\s*err_code\s*', '', target)
        elif pre_comma:
            print('Only pre')
            # err_code is last var in this definition block -> Keep indentation
            start, end = re.search(r',\s*err_code', target).span()
            try:
                new_l = re.sub(r',\s*err_code', ' '*(end-start), target)
            except:
                print('Pre-comma regular expression failed')
                sys.exit(1)
        elif post_comma:
            print('Only post')
            new_l = re.sub(r'err_code\s*,\s*', '', target)
        else:
            print('Regular expression extracting elements from line containing "err_code" failed: neither limiting character was found')
            sys.exit(1)
    except Exception as e:
        print('Regular expression extracting elements from line containing "err_code" failed: re.match returned nothing')
        sys.exit(1)
    return new_l

#end remove_var_in_line

#######################################################################################################################################
def replace_in_signature(lines):
    l_sig = False
    for i,l in enumerate(lines):
        if re.match(r'[^\!]*?(subroutine|function)\s*[a-Z][a-Z0-9_]*[\s]*[^\!][\(&]'):
            # Found start of signature
            if 'err_code' in l:
                l = re.sub('err_code', 'err_info', l)
                lines[i] = l
            else:
                l_sig = True

#end replace_in_signature

#######################################################################################################################################
def replace_in_directives(lines):
    # Replace err_code in acc directives
    # TODO: Include omp DONE
    out = []
    # Marker indicating that we are within an acc directive
    l_acc_omp = False
    # Return value indicating if anything was changed in this procedure
    l_modified = False
    for i, line in enumerate(lines):
        # If we are not already within a call, find the start of one
        if l_acc_omp == False:
            if line.strip().startswith('!$acc') or line.strip().startswith('!$omp'):
                print('ACC/OMP directive found in line ', i, ': ', line)
                combined_line, end_i = combine_lines(i, lines)
                print(combined_line)
                if 'err_code' in combined_line:
                    print('err_code in acc/omp directive found in line ', i, ': ', combined_line)
                    l_acc_omp = True
                    l_modified = True
        # Check for err_code starting on the same line
        if l_acc_omp == True and 'err_code' in line:
            print('err_code within acc/omp directive found in line ', i, ': ', line)
            l_acc_omp = False
            # Replace err_code with err_info
            out.append(line.replace('err_code', 'err_info'))
        else:
            # Otherwise, just append old line
            out.append(line)

    return out, l_modified

#end replace_in_directives

#######################################################################################################################################
def replace_in_checks_and_assignments(lines):
    # Replace err_code in checks and assignments, and add err_header to fstderr outputs
    out = []
    # Marker indicating that we are within an acc directive
    l_replace = False
    l_err_header = False
    # Return value indicating if anything was changed in this procedure
    l_modified = False
    for i, line in enumerate(lines):
        # If we have not found a check or assignment, find the start of one
        assignment_right = re.match(r'[!\s]*[^=]=\s*err_code', line)
        assignment_left = re.match(r'[!\s]*(?<!err_info%)err_code\s*=[^=]', line)
        check = re.search(r'(?<!err_info%)err_code\s*[/=]=', line)
        # TODO: Search for line to insert err_header output
        if check is not None or assignment_left is not None or assignment_right is not None:
            print('Check or assignment found in line ', i, ': ', line)
            combined_line, end_i = combine_lines(i, lines)
            print(combined_line)
            if 'err_code' in combined_line:
                print('err_code in check or assignment found in line ', i, ': ', combined_line)
                l_replace = True
                l_modified = True
        # Check for err_code starting on the same line
        # TODO: Do I need to make that 2nd check more specific?
        if l_replace == True and 'err_code' in line:
            print('err_code within check or assignment found in line ', i, ': ', line)
            l_replace = False
            # Replace err_code with err_info
            out.append(line.replace('err_code', 'err_info%err_code'))
        else:
            # Otherwise, just append old line
            out.append(line)
        # TODO: Fix 2nd check! Does not work -> If we can't find a following fstderr, then this whole loop collapses!
        if l_err_header and check is not None:
            if 'fstderr' in line:
                l_err_header = False
                # Get correct indentation from line
                indent = re.match(r'(\s*)[^\s]', line).groups()[0]
                out.append(indent + 'write(fstderr, *) err_header\n')
                out.append(line)

    return out, l_modified

#end replace_in_checks_and_assignments

#######################################################################################################################################
def modify_calls(lines):
    # TODO: Better differentiate the function search
    print('modify_calls')
    # Find all occurrences of err_code within a call to a function or a subroutine
    out = []
    # Marker indicating that we are within a procedure call
    l_call = False
    # Return value indicating if anything was changed in this procedure
    l_modified = False
    for i, line in enumerate(lines):
        # If we are not already within a call, find the start of one
        if l_call == False:
            if check_call(line):
                print('(Possible) Call found in line ', i, ': ', line)
                combined_line, end_i = combine_lines(i, lines)
                print(combined_line)
                # Check if line contains an unmodified 'err_code'
                if 'err_code' in combined_line and 'err_info%err_code' not in combined_line:
                    print('err_code in call found in line ', i, ': ', combined_line)
                    l_call = True
                    l_modified = True
        # Check for err_code starting on the same line
        if l_call == True and 'err_code' in line:
            print('err_code within call found in line ', i, ': ', line)
            l_call = False
            # Replace err_code with err_info
            out.append(line.replace('err_code', 'err_info'))
        else:
            # Otherwise, just append old line
            out.append(line)

    return out, l_modified

#end modify_calls

#######################################################################################################################################
def check_call(line):
    # Check of the given string represents the start of a call to a procedure
    # NOTE: We might not be able to distinguish between function calls and array slicing, but I don't think we care
    # 1. Check for subroutine call -> keyword call appears at the start of the line
    subroutine_check = re.match(r'[^!]*\bcall\b', line)
    # 2. Check for function call -> <func_name>
    # TODO: How do I decide between parentheses for calls and other parentheses????
    function_check = re.match(r'[^!]*\w+\s*\(', line)
    if subroutine_check is not None or function_check is not None:
        if subroutine_check is not None:
            print('Subroutine found')
            pass
        else:
            print('Function maybe found')
            pass

    return subroutine_check is not None or function_check is not None

#end check_call

#######################################################################################################################################
def find_use_blocks(fname):
    # NOT USED
    # List of all use blocks in file <fname>
    uses = {}
    rname = 'module'
    with open(fname, "r") as f:
        # Flag indicating if we are inside a use block
        # At the start of the file, we are not
        l_inside_use_block = False
        # Iterate through lines in file
        for l in f:
            print(l)
            # If inside use block
            if l.strip().startswith("subroutine") or l.strip().startswith("function"):
                rname = re.search(r'\s([^\s(]*)', l.strip()).groups()[0]
                print(rname)
            if l_inside_use_block:
                print("Inside block")
                # "implicit none" indicates end of use block
                if "implicit none" not in l:
                    print("No end")
                    # If not, just append line and continue
                    lines.append(l)
                else:
                    print("Block ended")
                    # Found end of block -> append entire block to uses
                    uses[rname] = "".join(lines)
                    # Empty lines list
                    lines = []
                    # Switch flag
                    l_inside_use_block = False
            else:
                print("Not inside use block")
                # If not inside use block, look out for use
                if l.strip().startswith("use"):
                    print("Use block start found")
                    l_inside_use_block = True
                    lines = [l]
        return uses

#######################################################################################################################################
if __name__ == "__main__":
    # parse the command line arguments
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('replace_paths', help="paths representing either a CLUBB source file or a directory containing such in any of the subfolders", nargs='*')

    args = parser.parse_args(sys.argv[1:])

    # call the main function
    main(args.replace_paths)
