"""Audit draft source contracts and expression structure; never execute physics."""
from pathlib import Path
import ast
import hashlib
import json
import re
import sys

folder = Path(__file__).resolve().parent
root = folder.parents[2]
sys.path.insert(0, str(root / 'clubb_kokkos/src/f2k'))
from fortran_procedures import find_procedure
from generate_derived_types import split_top

native = (root / 'src/CLUBB_core/saturation.F90').read_text()
target = (folder / 'saturation.hpp').read_text()
code = re.sub(r'//[^\n]*', '', target)
names = ['sat_mixrat_liq_k', 'sat_mixrat_liq_2D', 'sat_vapor_press_liq',
         'sat_vapor_press_liq_flatau', 'sat_vapor_press_liq_bolton',
         'sat_vapor_press_liq_gfdl', 'sat_vapor_press_liq_lookup',
         'sat_mixrat_ice_k', 'sat_mixrat_ice_2D', 'rcm_sat_adj']
blocks = {}
native_blocks = {}
for name in names:
    native_blocks[name] = re.search(
        r'^  (?:elemental )?(?:subroutine|function) ' + name +
        r'\(.*?^  end (?:subroutine|function) ' + name + r'\b', native, re.M | re.S).group()
    matches = list(re.finditer(r'\b' + name + r'\s*\(([^()]*)\)\s*\{', code))
    assert len(matches) == 1, (name, len(matches))
    match = matches[0]
    args = [arg.split('=')[0].strip().split()[-1].lstrip('&*')
            for arg in split_top(match[1])]
    assert args == [arg.name for arg in find_procedure(name).arguments], (name, args)
    depth, pos = 1, match.end()
    while depth:
        depth += (code[pos] == '{') - (code[pos] == '}')
        pos += 1
    blocks[name] = code[match.end():pos - 1]
assert [name for _, name in sorted((code.index(re.search(
    r'\b' + name + r'\s*\([^()]*\)\s*\{', code).group()), name) for name in names)] == names


def extra_block(name):
    match = re.search(r'\b' + name + r'\s*\([^()]*\)\s*\{', code)
    depth, pos = 1, match.end()
    while depth:
        depth += (code[pos] == '{') - (code[pos] == '}')
        pos += 1
    return code[match.end():pos - 1]


def logical_lines(text):
    pending = ''
    for line in text.splitlines():
        line = line.split('!')[0].strip()
        if not line:
            continue
        pending += ' ' + line.strip('& ')
        if not line.endswith('&'):
            yield pending.strip()
            pending = ''
    assert not pending


def assignments(text, variable, cpp=False):
    if cpp:
        return re.findall(r'\b' + variable + r'(?:\(i,\s*k\))?\s*=\s*([^;]+);', text)
    return [match[1] for line in logical_lines(text)
            if (match := re.match(variable + r'(?:\(i,k\))?\s*=\s*(.*)', line, re.I))]


def tree(expression, cpp=False):
    expression = re.sub(r'\s+', ' ', expression).strip()
    expression = expression.replace('Kokkos::', '').replace('saturation_detail::', '')
    expression = re.sub(r'\b(?:p_in_Pa|T_in_K|esat|esat_ice)\(i,k\)',
                        lambda m: m[0].split('(')[0], expression)
    expression = re.sub(r'pressure\(i,\s*k\)', 'p_in_Pa', expression)
    expression = expression.replace('excess > zero_threshold ? excess : zero_threshold',
                                    'max(excess, zero_threshold)')
    expression = re.sub(r'_core_rknd\b', '', expression)
    # Native kind-bearing literals and explicit C++ core casts have the same
    # value conversion boundary; this compares operation trees, not FP results.
    node = ast.parse(expression, mode='eval').body

    def norm(value):
        if isinstance(value, ast.Name):
            if value.id == 'excess':
                return tree('rtm - sat_mixrat_liq_api(p_in_Pa,theta*exner,saturation_formula)')
            return ('name', value.id.lower())
        if isinstance(value, ast.Constant):
            return ('constant', float(value.value).hex())
        if isinstance(value, ast.UnaryOp):
            return (type(value.op).__name__, norm(value.operand))
        if isinstance(value, ast.Call):
            fname = value.func.id.lower()
            if fname == 'real_t':
                return norm(value.args[0])
            if fname == 'pow':
                return ('power', norm(value.args[0]), norm(value.args[1]))
            return ('call', fname, tuple(norm(arg) for arg in value.args))
        if isinstance(value, ast.BinOp):
            if isinstance(value.op, ast.Pow):
                if isinstance(value.right, ast.Constant) and value.right.value == 2:
                    return ('Mult', norm(value.left), norm(value.left))
                return ('power', norm(value.left), norm(value.right))
            return (type(value.op).__name__, norm(value.left), norm(value.right))
        raise TypeError(ast.dump(value))
    return norm(node)


compared = {}

def compare(native_name, target_block, variables, coefficients=False):
    count = 0
    for variable in variables:
        original = assignments(native_blocks[native_name], variable)
        translated = assignments(target_block, variable, cpp=True)
        if coefficients:
            table = native_blocks[native_name]
            start = table.index('100._core_rknd * (/')
            values = re.findall(r'([\d.]+(?:[eE][+-]?\d+)?)_core_rknd',
                                table[start:table.index('/)', start)])[1:]
            original = [re.sub(r'a\((\d)\)',
                               lambda m: '(100._core_rknd * ' + values[int(m[1]) - 1] +
                               '_core_rknd)', expression) for expression in original]
        assert len(original) == len(translated), (native_name, variable, original, translated)
        for expected, actual in zip(original, translated):
            assert tree(expected) == tree(actual, True), (native_name, variable, expected, actual)
            count += 1
    compared[native_name] = count


compare('sat_mixrat_liq_k', extra_block('liquid_esat_for_mixrat'), ['T_in_C_sqd', 'esat'])
for name in ['sat_vapor_press_liq_flatau', 'sat_vapor_press_liq_bolton',
             'sat_vapor_press_liq_gfdl', 'sat_vapor_press_liq_lookup']:
    compare(name, blocks[name], ['T_in_C_sqd', 'esat'] if name.endswith('flatau') else ['esat'])
compare('sat_mixrat_ice_2D', blocks['sat_mixrat_ice_k'], ['esat_ice'], coefficients=True)
mixing_ratio_counts = {}
for source_name, target_name in [('sat_mixrat_liq_k', 'sat_mixrat_liq_k'),
                                  ('sat_mixrat_liq_2D', 'sat_mixrat_liq_2D'),
                                  ('sat_mixrat_ice_2D', 'sat_mixrat_ice_k')]:
    original = assignments(native_blocks[source_name], source_name)
    translated = assignments(blocks[target_name], 'result', cpp=True)
    assert len(original) == len(translated) == 3, (source_name, original, translated)
    for expected, actual in zip(original, translated):
        assert tree(expected) == tree(actual, True), (source_name, expected, actual)
    mixing_ratio_counts[source_name] = 3
compare('rcm_sat_adj', blocks['rcm_sat_adj'], ['answer', 'rcm'])

# The extracted scalar/2D liquid vapor-pressure blocks differ only by loops and
# array indexing. Compare each native assignment tree before sharing the body.
for variable in ['T_in_C_sqd', 'esat']:
    scalar = assignments(native_blocks['sat_mixrat_liq_k'], variable)
    array = assignments(native_blocks['sat_mixrat_liq_2D'], variable)
    assert len(scalar) == len(array)
    for expected, actual in zip(scalar, array):
        assert tree(expected) == tree(actual)

rawtable = native[native.index('    (/ 0.049560547'):native.index(' /)', native.index('    (/ 0.049560547'))]
values = re.findall(r'([\d.]+(?:[eE][+-]?\d+)?)_core_rknd', rawtable)
cpp_values = re.findall(r'case (\d+):\s*return real_t\(([^()]*)\);', extra_block('svp_liq_lookup_table'))
assert [int(index) for index, _ in cpp_values] == list(range(188, 344))
assert [value for _, value in cpp_values] == values
assert len(values) == 156

# Native public scalar SVP branch/call sequence remains explicit.
def call_names(text):
    return re.findall(r'\bsat_vapor_press_liq_(bolton|flatau|gfdl|lookup)\s*\(', text)
assert call_names(blocks['sat_vapor_press_liq']) == ['bolton', 'flatau', 'gfdl', 'lookup']
assert re.findall(r'case\s+(saturation_\w+)\s*:', extra_block('liquid_esat_for_mixrat')) == [
    'saturation_flatau', 'saturation_bolton', 'saturation_gfdl', 'saturation_lookup']
assert re.findall(r'case\s+(saturation_\w+)\s*:', blocks['sat_mixrat_ice_k']) == [
    'saturation_bolton', 'saturation_flatau', 'saturation_gfdl']
assert code.count('std::optional<int> start_index_in = std::nullopt') == 2
assert '*start_index_in - 1' in blocks['sat_mixrat_liq_2D']
assert 'stop_index = *start_index_in' in blocks['sat_mixrat_liq_2D']
assert code.count('std::is_same_v<typename Pressure::memory_space, typename Real2D::memory_space>') == 2
assert code.count('Kokkos::round(T_in_K)') == 2
assert code.count('Kokkos::WithoutInitializing') == 4
assert code.count('Kokkos::parallel_for(') == 2
assert 'Kokkos::abort("Error in rcm_sat_adj: could not determine rcm")' in code
assert 'if (iteration == itermax)' in blocks['rcm_sat_adj']
assert 'iteration <= itermax' in blocks['rcm_sat_adj']
assert 'f2k::' not in code
assert 'fma(' not in code
for left, right in [('(', ')'), ('[', ']'), ('{', '}')]:
    clean = re.sub(r'"(?:\\.|[^"\\])*"', '', code)
    depth = 0
    for char in clean:
        depth += (char == left) - (char == right)
        assert depth >= 0, (left, right)
    assert depth == 0, (left, right, depth)
assert all(len(line) <= 100 and line.rstrip() == line and '\t' not in line
           for line in target.splitlines())
record = {
    'status': 'Source-audited draft; not compiled or runtime validated',
    'authority': 'src/CLUBB_core/saturation.F90; all 1185 lines read',
    'native_sha256': hashlib.sha256(native.encode()).hexdigest(),
    'routine_order': names,
    'native_argument_counts': {name: len(find_procedure(name).arguments) for name in names},
    'all_native_signature_names_and_order': True,
    'formula_operation_trees_compared': compared,
    'mixing_ratio_operation_trees_compared': mixing_ratio_counts,
    'shared_liquid_scalar_array_assignment_trees': True,
    'lookup_literal_tokens_exact': 156,
    'selector_and_public_svp_call_order': True,
    'optional_start_index': 'Native 1-based, ascending/descending/absent and empty range retained',
    'output_memory': 'Pressure memory space via array_types owner',
    'scalar_device_apis': ['sat_mixrat_liq_api', 'sat_mixrat_ice', 'sat_vapor_press_liq', 'rcm_sat_adj'],
    'manual_audit': [
        'All 10 Native routine signatures, scientific descriptions/references and source order',
        'Liquid and ice pressure cap and both mixing-ratio/specific-humidity branches',
        'Distinct source liquid scalar/2D specific-humidity denominators',
        'Distinct GFDL liquid mixrat 173.15K and public SVP 203.15K floors',
        'Flatau source factored/Horner associations, coefficient scale and minimum temperatures',
        'All selectors and original sentinel defaults, including ice lookup sentinel',
        'ANINT tie-away rounding, exact lookup values and 188:343 clamping',
        'rcm_sat_adj tolerance/limit, update/exit/bracket order and Native loop-counter edge',
        'No stats, persistent mutable state, f2k lower calls or speculative FMA fixes',
        'Necessary scalar-device/array-fusion/table adaptations classified in INTEGRATION.md'
    ],
    'pending_shared_dependencies': ['model_flags.hpp: saturation_bolton/gfdl/flatau/lookup and I_sat_sphum'],
    'runtime': 'No compilation, unit tests, oracle tests or whole-model runs; root owns integration',
    'target_sha256': hashlib.sha256(target.encode()).hexdigest()
}
(folder / 'source_audit.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps({key: record[key] for key in ['status', 'routine_order',
      'formula_operation_trees_compared', 'lookup_literal_tokens_exact', 'target_sha256']}, indent=2))
