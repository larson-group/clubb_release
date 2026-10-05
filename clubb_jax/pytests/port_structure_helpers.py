"""Routine inventories shared by the Microphysics and SILHS source-contract checks."""

import ast
import re


def fortran_routines(source):
    # Join continuation lines before parsing declarations and internal procedures.
    code = re.sub(r'!.*', '', source.read_text())
    code = re.sub(r'&\s*\n\s*&?', ' ', code)
    pattern = r'^\s*(?:(?:recursive|pure|elemental|real(?:\([^)]*\))?|integer|type\([^)]*\))\s+)*(?:subroutine|function)\s+(\w+)\s*\((.*?)\)'
    routines = {}
    stack = []
    for line in code.splitlines():
        if re.match(r'\s*end\s+(subroutine|function)\b', line, re.I):
            stack.pop()
            continue
        match = re.match(pattern, line, re.I)
        if match:
            name, args = match.groups()
            stack.append(name.lower())
            routines[tuple(stack)] = (list(filter(None, map(str.strip, args.lower().split(',')))), {})
        elif stack and '::' in line:
            match = re.search(r'intent\s*\(\s*(inout|in|out)\s*\)', line, re.I)
            if match:
                declaration = re.sub(r'\([^()]*\)', '', line.split('::')[1])
                for name in declaration.lower().split(','):
                    routines[tuple(stack)][1][name.strip().split('=')[0].strip()] = match[1].lower()
    return routines


def python_routines(tree, prefix=()):
    result = {}
    for node in tree.body:
        if isinstance(node, ast.FunctionDef):
            key = prefix + (node.name.lower(),)
            result[key] = node
            result.update(python_routines(node, key))
    return result
