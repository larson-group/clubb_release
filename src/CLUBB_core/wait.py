import re
import shutil
from pathlib import Path

class Directive:
    """Represents an OpenACC directive."""
    
    def __init__(self, lines):
        self.raw_lines = lines
        self.acc_lines = self._extract_acc_lines(lines)
        self.leading_spaces = self._get_leading_spaces(self.acc_lines[0]) if self.acc_lines else ''
        self.acc_spacing = self._get_acc_spacing(self.acc_lines[0]) if self.acc_lines else ''
        self.clauses = self._parse_clauses()
    
    def _extract_acc_lines(self, lines):
        """Extracts only the lines containing !$acc from the raw lines."""
        return [line for line in lines if re.match(r'^\s*!\$acc', line, re.IGNORECASE)]
    
    def _get_leading_spaces(self, line):
        """Extracts leading spaces from the first !$acc line."""
        match = re.match(r'^(\s*)!\$acc', line, re.IGNORECASE)
        return match.group(1) if match else ''
    
    def _get_acc_spacing(self, line):
        """Extracts the spaces immediately after !$acc."""
        match = re.match(r'^\s*!\$acc(\s+)', line, re.IGNORECASE)
        return match.group(1) if match else ''
    
    def _parse_clauses(self):
        """Parses the combined clauses of all !$acc lines."""
        combined_line = ''.join(line.strip() for line in self.acc_lines).rstrip('&')
        match = re.match(r'^\s*!\$acc\s+(.*)', combined_line, re.IGNORECASE)
        if match:
            return match.group(1).split()
        return []
    
    def contains_clause(self, clause):
        """Checks if the directive contains a specific clause."""
        return clause in self.clauses
    
    def contains_keyword(self, keyword):
        """Checks if the directive contains a specific keyword."""
        return any(keyword in clause for clause in self.clauses)
    
    def add_clause(self, clause):
        """Adds a clause to the directive if not already present."""
        if not self.contains_clause(clause):
            self.clauses.append(clause)
    
    def reconstruct(self):
        """Reconstructs the directive, appending wait to the last line."""
        if not self.acc_lines:
            return self.raw_lines  # No !$acc lines, return raw unmodified

        # Add wait to the last line of the directive
        last_acc_line = self.raw_lines[-1].rstrip('\n')
        if last_acc_line.endswith('&'):
            last_acc_line = last_acc_line.rstrip('&') + " wait &\n"
        else:
            last_acc_line = last_acc_line + " wait\n"

        # Replace the last line in the raw lines with the modified version
        modified_lines = self.raw_lines[:]
        modified_lines[-1] = last_acc_line
        return modified_lines

def read_file(file_path):
    with open(file_path, 'r') as file:
        return file.readlines()

def write_file(file_path, lines):
    with open(file_path, 'w') as file:
        file.writelines(lines)

def backup_file(original_file):
    backup_path = original_file.with_suffix('.bak')
    shutil.move(original_file, backup_path)
    return backup_path

def find_directives(lines):
    """Identifies OpenACC directives and their spans in the file."""
    directive_spans = []
    start_idx = None

    for idx, line in enumerate(lines):
        if re.match(r'^\s*!\$acc', line, re.IGNORECASE):
            if start_idx is None:
                start_idx = idx
        if start_idx is not None and not line.strip().endswith('&'):
            # Directive ends when a line does not end with '&'
            directive_spans.append((start_idx, idx + 1))
            start_idx = None

    if start_idx is not None:  # Catch the last directive if unclosed
        directive_spans.append((start_idx, len(lines)))

    return directive_spans

def process_directive(lines):
    """Adds wait to every directive."""
    directive = Directive(lines)

    # Skip directives with async, end, update, or routine
    if (
        directive.contains_keyword("end") or 
        directive.contains_keyword("declare") or 
        directive.contains_keyword("routine")
    ):
        return lines

    directive = Directive(lines)
    directive.add_clause("wait")
    return directive.reconstruct()

def process_file(input_file):
    """Processes the input file to modify directives."""
    original_lines = read_file(input_file)
    #backup_file(input_file)
    modified_lines = original_lines.copy()

    directive_spans = find_directives(original_lines)
    for start, end in directive_spans:
        original_directive = original_lines[start:end]
        modified_directive = process_directive(original_directive)
        modified_lines[start:end] = modified_directive

    write_file(input_file, modified_lines)
    print(f"Processed file saved as {input_file}. Original backed up as {input_file}.bak")

# Entry point
if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python script.py <fortran_file>")
        sys.exit(1)

    input_file = Path(sys.argv[1])
    if not input_file.exists():
        print(f"Error: File {input_file} does not exist.")
        sys.exit(1)

    process_file(input_file)
