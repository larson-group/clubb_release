import re
import sys
import argparse

# Horizontal line separator (for output formatting)
horizontal_line = "-" * 100


def strip_strings(line: str) -> str:
    """Remove Fortran string literals delimited by ' or "."""
    out = []
    in_single = in_double = False
    for ch in line:
        if in_single:
            if ch == "'":
                in_single = False
            continue
        if in_double:
            if ch == '"':
                in_double = False
            continue
        if ch == "'":
            in_single = True
            continue
        if ch == '"':
            in_double = True
            continue
        out.append(ch)
    return "".join(out)

def check_implicit_none(filename, lines, verbose=False):
#
#     Description: This verifies that the number of program 
#     block and implicit none statements in a FORTRAN file match. 
#
#        This subroutine works by testing whether
#           ( number of programs,modules,subroutines, and functions )  
#                     ==  ( number of 'implicit none' statements ).
#
#     Returns
#     	True if the file has no missing implicit none statements.
#####################################################################

    current_stmt_line = None
    waiting_for_implicit = False
    warnings = []

    for linenum, line in enumerate(lines, start=1):
        stripped = line.strip().lower()

        # Detect start of a block
        if (stripped.startswith("function ")
            or stripped.startswith("subroutine ")
            or (stripped.startswith("module ") and not stripped.startswith("module procedure"))
            or stripped.startswith("program ")):
            
            if waiting_for_implicit:
                # We found a new block but never saw an implicit none for the last one
                warnings.append(f"line {current_stmt_line[0]}: {current_stmt_line[1].rstrip()}")
            
            # Start tracking this new statement
            current_stmt_line = (linenum, line)
            waiting_for_implicit = True
            if verbose:
                print(f"Found block start at line {linenum}: {line.rstrip()}")

        # Detect implicit none
        elif stripped.startswith("implicit none"):
            if waiting_for_implicit:
                if verbose:
                    print(f"Found implicit none for block starting at line {current_stmt_line[0]}")
                waiting_for_implicit = False
                current_stmt_line = None

    # If file ended and we were still waiting
    if waiting_for_implicit and current_stmt_line:
        warnings.append(f"line {current_stmt_line[0]}: {current_stmt_line[1].rstrip()}")

    # Reporting
    if warnings:
        print(horizontal_line)
        print(f"--------------------- Implicit None Check: {filename} ---------------------")
        print("WARNING: Missing 'implicit none' after these blocks:")
        for w in warnings:
            print(w)
        return False
    return True

def check_use_only(filename, lines, verbose=False):
#
#     Description: This verifies that the file does not contain
#     any "use" statements without being followed by an "only" restriction.
#
#     	Basically, this checks for any line of the form
#       	use module_name  ! Comment
#     	and if so, returns an error because there is no "only" statement.
#     
#     Returns
#     	True if the file has no unrestricted "use" statements.
####################################################################

    use_regex           = re.compile(r"^ *use +\w+ *(?:!.*)?$", re.IGNORECASE)
    offending_lines = []
    
    for linenum, line in enumerate(lines, start=1):
        if use_regex.match(line):
            offending_lines.append((linenum, line))

    if offending_lines:

        print(f"--------------------- Implicit None Check: {filename} ---------------------")
        print("WARNING: 'use' statements without 'only' found in the following lines:")
        print("WARNING: Check that comma is on same line as 'use', as CLUBB requires.")

        for linenum, content in offending_lines:
            # Print line number and content (content already has its newline)
            print(f"{linenum} : {content}", end="")

        return False
    return True

def check_private_default(filename, lines, verbose=False):
#
#     Description: This verifies that a file's modules
#     have a corresponding private statement to set the default
#     scope of the module to private.
#
#     This subroutine works by testing whether
#     	( number of module statements )  ==  ( number of private statements ).
#
#     Returns
#     	True if the module is set to private default scope.  
####################################################################

    current_module = None
    waiting_for_private = False
    warnings = []

    for linenum, line in enumerate(lines, start=1):
        stripped = line.strip().lower()

        # Detect module start (ignore "module procedure")
        if stripped.startswith("module ") and not stripped.startswith("module procedure"):
            if waiting_for_private:
                warnings.append(
                    f"line {current_module[0]}: {current_module[1].rstrip()}"
                )
            current_module = (linenum, line)
            waiting_for_private = True
            if verbose:
                print(f"Found module at line {linenum}: {line.rstrip()}")

        # Detect default private
        elif waiting_for_private:
            # Allow "private" or "private   ! comment", but reject "private :: ..."
            tokens = stripped.split()
            if tokens and tokens[0] == "private":
                if len(tokens) == 1 or tokens[1].startswith("!"):
                    if verbose:
                        print(f"Found default private for module starting at line {current_module[0]}")
                    waiting_for_private = False
                    current_module = None

    # End-of-file check
    if waiting_for_private and current_module:
        warnings.append(
            f"line {current_module[0]}: {current_module[1].rstrip()}"
        )

    if warnings:
        print(f"--------------------- Default Private Check: {filename} ---------------------")
        print("WARNING: Missing default 'private' in the following modules:")
        for w in warnings:
            print(w)
        return False
    return True

def check_line_length(filename, lines):
#
#     Description: This subroutine  verifies that no line in the file
#     is longer than the maxLength. 
#
#        This subroutine works by comparing the length of each line
#        to the maxLength
#
#     Returns
#     	True if no line in the file is greater than the maxLength.
#####################################################################

    max_length=100
    long_lines = []

    for linenum, line in enumerate(lines, start=1):
        length = len(line.rstrip("\n"))
        if length > max_length:
            long_lines.append((linenum, length, line))

    if long_lines:

        print(f"--------------------- Line Length Check: {filename} ---------------------")
        print(f" WARNING: Lines exceed {max_length} characters (limit {max_length}):")

        for linenum, length, content in long_lines:
            print(f"line {linenum} : {length} chars : {content.rstrip()}")

        return False

    return True


def check_forbidden_tokens(filename, lines, verbose=False):
#
#     Description: This subroutine verifies that there are no forbidden
#       tokens used in the Fortran source file.
#
#     Returns
#     	True if no line has forbidden tokens.
#####################################################################

    forbidden_regex     = re.compile(r"\.le\.|\.ge\.|\.lt\.|\.gt\.|\.eq\.|\.ne\.", re.IGNORECASE)
    issues = []

    for linenum, line in enumerate(lines, start=1):
        # Ignore any text after a '!' (comment) for this check
        code_part = line.split("!")[0]
        matches = list(forbidden_regex.finditer(code_part))
        if matches:
            issues.append((linenum, matches, line))

    if issues:
        print(f"--------------------- Check Forbidden Tokens: {filename} ---------------------")
        print("WARNING: Forbidden tokens found in the following lines:")

        for linenum, matches, content in issues:
            highlighted = content
            # apply highlights (from end to start so indices don’t shift)
            for m in reversed(matches):
                start, end = m.span()
                token = content[start:end]
                highlighted = (
                    highlighted[:start]
                    + f"\033[91m{token}\033[0m"
                    + highlighted[end:]
                )
            token_list = ", ".join({m.group(0).lower() for m in matches})
            print(f"line {linenum} : contains {token_list} : {highlighted.strip()} ")

        return False

    return True

def check_end_block(filename, lines, verbose=True):
#
#     Description: This subroutine verifies that the name in an end statment matches the
#     		name of the block.
#    
#              ( e.g.
#              		subroutine do_calc()
#              		end subroutine do_calc ! Passed!
#
#              		subroutine do_calc()
#              		end subroutine         ! Failed!
#              	)
#
#     Returns
#     	True if all names at the beginning of a block match their end statment.
#####################################################################

    # Tiny regex to extract a valid Fortran identifier at the start of a token
    identifier = re.compile(r"^[A-Za-z_]\w*")

    stack = []   # holds tuples of (type, name, line_number)
    warnings = []
    valid = True

    if verbose:
        print("----------------End Block Test----------------")

    for idx, rawline in enumerate(lines, start=1):
        # Strip inline comments
        code = rawline.split("!", 1)[0]
        code = strip_strings(code).strip()
        if not code:
            continue
        lower_line = code.lower()

        # ---- Handle END ----
        if lower_line.startswith("end"):
            parts = lower_line.split()
            if len(parts) >= 2 and parts[1] in ("program", "module", "function", "subroutine"):
                end_type = parts[1]
                end_name = parts[2] if len(parts) >= 3 else ""
                # normalize name if present
                if end_name:
                    m = identifier.match(end_name)
                    if m:
                        end_name = m.group(0)

                if not stack:
                    warnings.append(f"Line {idx}: Found 'end {end_type}' with no matching block.")
                    valid = False
                    continue

                start_type, start_name, start_line = stack.pop()

                if end_type != start_type:
                    warnings.append(
                        f"Line {idx}: Found 'end {end_type} {end_name}' but top of stack is "
                        f"'{start_type} {start_name}' (opened at line {start_line})."
                    )
                    valid = False
                elif not end_name:
                    warnings.append(
                        f"Line {idx}: 'end {end_type}' missing name, expected '{start_name}' "
                        f"(opened at line {start_line})."
                    )
                    valid = False
                elif end_name.lower() != start_name.lower():
                    warnings.append(
                        f"Line {idx}: 'end {end_type} {end_name}' does not match "
                        f"'{start_type} {start_name}' (opened at line {start_line})."
                    )
                    valid = False
                else:
                    if verbose:
                        print(f"Found matching end for {end_name}")
            continue

        # ---- Handle FUNCTION ----
        if "function" in lower_line:
            parts = lower_line.split()
            if "function" in parts:
                i = parts.index("function")
                if i + 1 < len(parts):
                    match = identifier.match(parts[i + 1])
                    if match:
                        name = match.group(0)
                        stack.append(("function", name, idx))
                        if verbose:
                            print(f"Found Function {name} at line {idx}")
                    continue

        # ---- Handle SUBROUTINE ----
        if "subroutine" in lower_line:
            parts = lower_line.split()
            if "subroutine" in parts:
                i = parts.index("subroutine")
                if i + 1 < len(parts):
                    match = identifier.match(parts[i + 1])
                    if match:
                        name = match.group(0)
                        stack.append(("subroutine", name, idx))
                        if verbose:
                            print(f"Found Subroutine {name} at line {idx}")
                    continue

        # ---- Handle MODULE ----
        if lower_line.startswith("module "):
            parts = lower_line.split()
            if len(parts) >= 2 and parts[1] != "procedure":
                match = identifier.match(parts[1])
                if match:
                    name = match.group(0)
                    stack.append(("module", name, idx))
                    if verbose:
                        print(f"Found Module {name} at line {idx}")
                continue

        # ---- Handle PROGRAM ----
        if lower_line.startswith("program "):
            parts = lower_line.split()
            if len(parts) >= 2:
                match = identifier.match(parts[1])
                if match:
                    name = match.group(0)
                    stack.append(("program", name, idx))
                    if verbose:
                        print(f"Found Program {name} at line {idx}")
                continue

    if verbose and not stack:
        print("stack empty -- no unmatched ends")

    # If stack not empty at the end → unmatched opens
    while stack:
        start_type, start_name, start_line = stack.pop()
        warnings.append(f"Unclosed '{start_type} {start_name}' starting at line {start_line}.")
        valid = False

    if not valid:
        print(f"--------------------- Check End Blocks: {filename} ---------------------")
        print("WARNING: End block mismatches detected")
        for w in warnings:
            print(w)
        return False

    return True

def main():
    parser = argparse.ArgumentParser(description="CLUBB Standards Check (Python version)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose messages during checks")
    parser.add_argument("files", nargs="+", help="Fortran source file(s) to check")
    args = parser.parse_args()
    verbose = args.verbose
    fail_count = 0

    print("CLUBBStandardsCheck.py has begun.")

    for filename in args.files:

        try:
            with open(filename, 'r') as f:
                lines = f.readlines()
        except IOError:
            print(f"Bad Filename: {filename}")
            fail_count += 1
            continue

        # # Run all checks on this file:
        if not check_implicit_none(filename, lines):
            print(horizontal_line)
            fail_count += 1

        if not check_use_only(filename, lines, verbose):
            print(horizontal_line)
            fail_count += 1

        if not check_private_default(filename, lines, verbose):
            print(horizontal_line)
            fail_count += 1

        if not check_forbidden_tokens(filename, lines, verbose):
            print(horizontal_line)
            fail_count += 1

        if not check_line_length(filename, lines):
            print(horizontal_line)
            fail_count += 1

        if not check_end_block(filename, lines, verbose):
            print(horizontal_line)
            fail_count += 1

    print("CLUBBStandardsCheck.py has finished.")
    if fail_count > 0:
        print(f"FAIL: {fail_count} check(s) failed.")
        sys.exit(1)
    else:
        print("PASS: 0 checks failed")
        sys.exit(0)

if __name__ == "__main__":
    main()
