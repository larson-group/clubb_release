import re
import sys
import argparse
import os

# Horizontal line separator (for output formatting)
horizontal_line = "-" * 100
CLUBB_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
identifier_regex = re.compile(r"^[A-Za-z_]\w*$")


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

def strip_comments_and_strings(line: str) -> str:
    """Remove string literals and trailing Fortran comments from a line."""
    return strip_strings(line).split("!", 1)[0].rstrip("\n")

def iter_logical_statements(lines):
    """Yield free-form Fortran logical statements as (start_line, end_line, stmt)."""
    buffer = []
    start_line = None

    for linenum, raw_line in enumerate(lines, start=1):
        code = strip_comments_and_strings(raw_line)
        if not code.strip():
            continue

        segment = code.rstrip()
        continued = segment.endswith("&")
        if continued:
            segment = segment[:-1].rstrip()

        if not buffer:
            start_line = linenum
        else:
            segment = segment.lstrip()
            if segment.startswith("&"):
                segment = segment[1:].lstrip()

        if segment:
            buffer.append(segment)

        if not continued:
            yield start_line, linenum, " ".join(buffer).strip()
            buffer = []
            start_line = None

    if buffer:
        yield start_line, start_line, " ".join(buffer).strip()

def parse_use_only_imports(stmt):
    """Return imported local names for a `use ..., only:` statement."""
    lower_stmt = stmt.lower()
    if not lower_stmt.startswith("use"):
        return []
    if " only:" not in lower_stmt:
        return []

    only_part = stmt.split("only:", 1)[1]
    imports = []
    for item in only_part.split(","):
        name = item.strip()
        if not name:
            continue
        if "=>" in name:
            name = name.split("=>", 1)[0].strip()
        lower_name = name.lower()
        if identifier_regex.match(lower_name):
            imports.append(lower_name)
    return imports

def check_threadprivate_vars(filename, lines, verbose=False):
    #
    #     Description: This subroutine verifies that module variables in
    #       CLUBB_core and SILHS source files are declared threadprivate,
    #       and that threadprivate directives do not reference undeclared
    #       module variables.
    #
    #     Returns
    #       True if the file's module declarations and threadprivate lists match.
    #####################################################################

    top_level_block_type = None
    for _, _, stmt in iter_logical_statements(lines):
        lower_stmt = stmt.lower()
        if lower_stmt.startswith("module ") and not lower_stmt.startswith("module procedure"):
            top_level_block_type = "module"
            break
        if lower_stmt.startswith("program "):
            top_level_block_type = "program"
            break

    if top_level_block_type != "module":
        return True

    def iter_declaration_statements():
        parts = []
        in_declaration = False
        in_type_declaration = False

        for raw_line in lines:
            line = strip_comments_and_strings(raw_line)
            if not line.strip():
                continue

            if re.search(r"\bcontains\b", line, re.IGNORECASE):
                break

            if (
                re.search(r"\btype\b", line, re.IGNORECASE) is not None
                and re.search(r"\bend\s*\btype\b", line, re.IGNORECASE) is None
                and re.search(r"\btype\s*\(", line, re.IGNORECASE) is None
            ):
                in_type_declaration = True
                continue

            if re.search(r"\bend\s*\btype\b", line, re.IGNORECASE) is not None:
                in_type_declaration = False
                continue

            is_decl_start = any(
                re.search(pattern, line, re.IGNORECASE) is not None
                for pattern in (
                    r"\breal\b",
                    r"\binteger\b",
                    r"\blogical\b",
                    r"\bcharacter\b",
                    r"\btype\s*\(",
                )
            )

            if (not is_decl_start and not in_declaration) or in_type_declaration:
                continue

            if "&" in line:
                split_line = line.split("&")
                parts.append(split_line[0])
                in_declaration = True
            else:
                yield "".join(parts) + line
                parts = []
                in_declaration = False

        if parts:
            yield "".join(parts)

    def find_declared_module_vars():
        declared = []
        blacklist = {"reshape", "kind", "dimension", "public", "private"}

        def split_declarators(decl):
            items = []
            current = []
            depth = 0
            for ch in decl:
                if ch == "(":
                    depth += 1
                elif ch == ")":
                    depth = max(depth - 1, 0)
                elif ch == "," and depth == 0:
                    item = "".join(current).strip()
                    if item:
                        items.append(item)
                    current = []
                    continue
                current.append(ch)
            item = "".join(current).strip()
            if item:
                items.append(item)
            return items

        for declaration_stmt in iter_declaration_statements():
            if re.search(r"\bparameter\b(?=.*::)", declaration_stmt, re.IGNORECASE):
                continue
            if "::" not in declaration_stmt:
                continue

            decl = declaration_stmt.split("::", 1)[1]
            candidates = []
            for item in split_declarators(decl):
                item = item.split("=", 1)[0].strip()
                match = re.match(r"([A-Za-z_]\w*)", item)
                if match:
                    candidates.append(match.group(1))
            filtered = [name for name in candidates if name.lower() not in blacklist]
            filtered = [name for name in filtered if not name.endswith("_c")]
            filtered = [name for name in filtered if not name.startswith("_")]
            filtered = [name for name in filtered if name not in {"d_var_total"}]
            declared.extend(filtered)

        return declared

    def find_threadprivate_vars():
        src_string = "".join(lines)
        src_string = re.split(r"contains\s*\n", src_string, flags=re.IGNORECASE)[0]
        threadprivate_groups = re.findall(
            r"(?si)!\$omp\s*threadprivate\s*(\(.*?\))",
            src_string,
        )
        flattened = []
        for group in threadprivate_groups:
            for entry in group.strip("()").split(","):
                cleaned = re.sub(r"(?i)&\s*!\$omp", "", entry).strip()
                if cleaned:
                    flattened.append(cleaned)
        return flattened

    declared_vars = set(find_declared_module_vars())
    threadprivate_vars = set(find_threadprivate_vars())

    missing_threadprivate = sorted(declared_vars - threadprivate_vars)
    missing_declaration = sorted(threadprivate_vars - declared_vars)

    if missing_threadprivate or missing_declaration:
        print(f"--------------------- Threadprivate Check: {filename} ---------------------")
        if missing_threadprivate:
            print("WARNING: Missing threadprivate for:")
            for name in missing_threadprivate:
                print(f"  {name}")
        if missing_declaration:
            print("WARNING: Threadprivate without declaration for:")
            for name in missing_declaration:
                print(f"  {name}")
        return False

    return True

def check_unused_use_imports(filename, lines, verbose=False):
    #
    #     Description: This subroutine verifies that names imported via
    #       "use ..., only:" are referenced within their containing scope.
    #
    #     Module-level imports may be used anywhere in the module, including
    #       contained procedures. Procedure-level imports must be used within
    #       that procedure.
    #
    #     Returns
    #       True if all imported names are referenced.
    #####################################################################

    def push_scope(scope_type, name, line):
        scope = {
            "type": scope_type,
            "name": name.lower(),
            "line": line,
            "imports": [],
            "used_names": set(),
        }
        scope_stack.append(scope)
        return scope

    scope_stack = []
    warnings = []

    for start_line, end_line, stmt in iter_logical_statements(lines):
        lower_stmt = stmt.lower()

        if lower_stmt == "contains":
            continue

        if lower_stmt.startswith("end "):
            parts = lower_stmt.split()
            if len(parts) >= 2 and parts[1] in ("module", "program", "subroutine", "function"):
                if scope_stack:
                    scope = scope_stack.pop()
                    for imported in scope["imports"]:
                        if imported["local_name"] not in scope["used_names"]:
                            warnings.append(
                                (
                                    imported["line"],
                                    scope["type"],
                                    scope["name"],
                                    imported["local_name"],
                                    imported["stmt"],
                                )
                            )
            continue

        if lower_stmt.startswith("module ") and not lower_stmt.startswith("module procedure"):
            name = lower_stmt.split()[1]
            push_scope("module", name, start_line)
            continue

        if lower_stmt.startswith("program "):
            name = lower_stmt.split()[1]
            push_scope("program", name, start_line)
            continue

        if "subroutine " in lower_stmt:
            match = re.search(r"\bsubroutine\s+([a-z_]\w*)", lower_stmt)
            if match:
                push_scope("subroutine", match.group(1), start_line)
                continue

        if "function " in lower_stmt:
            match = re.search(r"\bfunction\s+([a-z_]\w*)", lower_stmt)
            if match:
                push_scope("function", match.group(1), start_line)
                continue

        imported_names = parse_use_only_imports(stmt)
        if imported_names:
            if scope_stack:
                target_scope = scope_stack[-1]
            else:
                target_scope = push_scope("file", os.path.basename(filename), start_line)
            for local_name in imported_names:
                target_scope["imports"].append(
                    {"local_name": local_name, "line": start_line, "stmt": stmt}
                )
            continue

        tokens = {
            token.lower()
            for token in re.findall(r"\b[a-z_]\w*\b", lower_stmt)
        }
        tokens.update(
            token.lower()
            for token in re.findall(r"(?<=_)([a-z_]\w*)\b", lower_stmt)
        )
        if not tokens:
            continue

        for scope in scope_stack:
            scope["used_names"].update(tokens)

    while scope_stack:
        scope = scope_stack.pop()
        for imported in scope["imports"]:
            if imported["local_name"] not in scope["used_names"]:
                warnings.append(
                    (
                        imported["line"],
                        scope["type"],
                        scope["name"],
                        imported["local_name"],
                        imported["stmt"],
                    )
                )

    if warnings:
        print(f"--------------------- Unused Use Imports Check: {filename} ---------------------")
        print("WARNING: Unused names imported by 'use ..., only:' found:")
        for line, scope_type, scope_name, local_name, stmt in warnings:
            print(
                f"line {line} : unused import '{local_name}' in {scope_type} "
                f"'{scope_name}' : {stmt}"
            )
        return False

    return True

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
    in_module_contains = False

    for linenum, line in enumerate(lines, start=1):
        stripped = line.strip().lower()

        if stripped == "contains":
            in_module_contains = True
            waiting_for_implicit = False
            current_stmt_line = None
            continue

        # Detect start of a block
        if (stripped.startswith("function ")
            or stripped.startswith("subroutine ")
            or (stripped.startswith("module ") and not stripped.startswith("module procedure"))
            or stripped.startswith("program ")):
            if in_module_contains and (
                stripped.startswith("function ") or stripped.startswith("subroutine ")
            ):
                continue
            
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
        stripped = line.lstrip()
        # Ignore comment-only lines and accelerator/OpenMP directives.
        if stripped.startswith("!") or stripped.startswith("#"):
            continue
        code_only = strip_strings(line).split("!", 1)[0].rstrip("\n").rstrip()
        length = len(code_only)
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

def format_target_label(files):
    """Return a concise label describing the current check target(s)."""
    parent_dirs = sorted({os.path.dirname(os.path.abspath(path)) for path in files})
    rel = lambda path: os.path.relpath(path, start=CLUBB_ROOT)

    if len(parent_dirs) == 1:
        return rel(parent_dirs[0])

    common_dir = os.path.commonpath(parent_dirs)
    if common_dir and common_dir != os.path.sep:
        return f"multiple directories under {rel(common_dir)}"

    return "multiple directories"

def main():
    parser = argparse.ArgumentParser(description="CLUBB Standards Check (Python version)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose messages during checks")
    parser.add_argument("files", nargs="+", help="Fortran source file(s) to check")
    args = parser.parse_args()
    verbose = args.verbose
    fail_count = 0
    target_label = format_target_label(args.files)

    print(f"CLUBBStandardsCheck.py: beginning check on {target_label}")

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

        if not check_unused_use_imports(filename, lines, verbose):
            print(horizontal_line)
            fail_count += 1

        if not check_threadprivate_vars(filename, lines, verbose):
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

    if fail_count > 0:
        print(f"FAIL: {fail_count} check(s) failed for {target_label}")
        sys.exit(1)
    else:
        print(f"PASS: 0 checks failed for {target_label}")
        sys.exit(0)

if __name__ == "__main__":
    main()
