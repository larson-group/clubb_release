import re
import shutil
from pathlib import Path

has_pending_wait = False

class Directive:
    """Represents an OpenACC directive."""
    
    def __init__(self, lines):
        self.raw_lines = lines
        self.acc_lines = self._extract_acc_lines(lines)
        self.clauses = self._parse_clauses_with_formatting()
        print(self.clauses)
    
    def _extract_acc_lines(self, lines):
        """Extracts only the lines containing !$acc from the raw lines."""
        return [line for line in lines if re.match(r'^\s*!\$acc', line, re.IGNORECASE)]
    
    def _parse_clauses_with_formatting(self):
        """Parses clauses along with spaces and line breaks."""
        clauses = []
        for line in self.acc_lines:
            leading_spaces = re.match(r'^(\s*)', line).group(1)  # Capture leading spaces
            clauses.append(leading_spaces)
            
            match = re.match(r'^\s*(!\$acc\s+)(.*)(&\s+$)?', line, re.IGNORECASE)
            if match:
                clauses.append(match.group(1))  # !$acc with spaces
                clauses.extend(match.group(2).split())  # Split the rest into clauses
                if match.lastindex >= 3:  
                    clauses.extend(match.group(3))  # Append the "&" explicitly or handle as needed
            
            clauses.append("\n")  # End of directive line

        return clauses
    
    def contains_clause(self, clause):
        """Checks if the directive contains a specific clause."""
        return any(c == clause for c in self.clauses)
    
    def contains_keyword(self, keyword):
        """Checks if the directive contains a specific keyword."""
        return any(keyword in clause for clause in self.clauses)
    
    def add_clause(self, clause):
      """Adds a clause to the directive just before the last newline."""
      if not self.contains_clause(clause):
          # Find the last occurrence of "\n" or "&\n" and insert before it
          for i in reversed(range(len(self.clauses))):
              if self.clauses[i] in ["\n"]:
                  self.clauses.insert(i, clause)
                  break
      print(self.clauses)
    
    def remove_clause(self, clause_prefix):
        """Removes any clause that starts with the given prefix."""
        self.clauses = [clause for clause in self.clauses if not clause.startswith(clause_prefix)]
    
    def reconstruct(self):
        """Reconstructs the directive using stored clauses."""
        reconstructed_lines = []
        current_line = ""
        for clause in self.clauses:
            if clause in ["\n"]:  # Line break or continuation
                current_line += clause
                reconstructed_lines.append(current_line)
                current_line = ""
            elif clause.strip() and "acc" not in clause:  # Only add a space if the clause is not all spaces
                current_line += clause + " "
            else:
                current_line += clause

        if current_line.strip():  # Add any remaining content
            reconstructed_lines.append(current_line + "\n")

        return reconstructed_lines

    def get(self, n):
        if 0 <= n <= len(self.clauses)-1:
            return self.clauses[n]
        return None

def read_file(file_path):
    with open(file_path, 'r') as file:
        return file.readlines()

def write_file(file_path, lines):
    with open(file_path, 'w') as file:
        file.writelines(lines)

# def backup_file(original_file):
#     backup_path = original_file.with_suffix('.bak')
#     shutil.move(original_file, backup_path)
#     return backup_path

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

def process_directive(lines, has_pending_wait=False):
    """Processes a directive, modifying it as needed."""
    directive = Directive(lines)

    # Skip directives with async, end, update, or routine
    if (
        has_pending_wait and 
        directive.contains_keyword("end") and directive.contains_keyword("parallel")
    ):
        directive.add_clause("\n"+directive.get(0)+"!$acc wait")
        has_pending_wait = False
    elif  (
        directive.contains_keyword("reduction")
    ):
        directive.remove_clause("wait")
        directive.add_clause("async(1)")
        has_pending_wait = True  # Mark that a `wait` needs to be added later
    elif  (
        directive.contains_keyword("copy") or 
        directive.contains_keyword("copyout") or 
        directive.contains_keyword("exit") or 
        directive.contains_keyword("update")
    ):
        directive.remove_clause("async(1)")
        directive.add_clause("wait")
    elif (
        directive.contains_keyword("end") or 
        directive.contains_keyword("declare") or 
        directive.contains_keyword("routine") or 
        directive.contains_keyword("async")
    ):
        return lines, has_pending_wait
    elif  (
        directive.contains_keyword("enter") or 
        directive.contains_keyword("data") or
        ( directive.contains_keyword("parallel") and directive.contains_keyword("loop") )
    ):
        directive.remove_clause("wait")
        directive.add_clause("async(1)")
        

    return directive.reconstruct(), has_pending_wait

def process_file(input_file):
    """Processes the input file to modify directives."""
    original_lines = read_file(input_file)
    #backup_file(input_file)
    modified_lines = original_lines.copy()

    has_pending_wait = False
    directive_spans = find_directives(original_lines)
    for start, end in directive_spans:
        original_directive = original_lines[start:end]
        modified_directive, has_pending_wait = process_directive(original_directive, has_pending_wait)
        print(f"has_pending_wait = {has_pending_wait}")
        modified_lines[start:end] = modified_directive

    write_file(input_file, modified_lines)
    print(f"Processed file saved as {input_file}.")# Original backed up as {input_file}.bak")

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
