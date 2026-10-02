"""Structural audit of the source-shaped microphysics interface (cores exempt)."""
import ast
import re
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[2]
MODULES = sorted((ROOT / 'src/Microphys').glob('*.F90'))


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


def test_interface_file_inventory():
    expected = {p.stem for p in MODULES}
    target = ROOT / 'clubb_jax/src/Microphys'
    assert {p.stem for p in target.glob('*.py') if p.stem != '__init__'} == expected
    assert {p.name for p in target.iterdir() if p.is_dir() and p.name != '__pycache__'} <= {
        p.name for p in (ROOT / 'src/Microphys').iterdir() if p.is_dir()}


@pytest.mark.parametrize('source', MODULES, ids=lambda p: p.stem)
def test_source_inventory_order_and_argument_order(source):
    target = ROOT / 'clubb_jax' / source.relative_to(ROOT).with_suffix('.py')
    assert target.exists()
    routines = fortran_routines(source)
    functions = python_routines(ast.parse(target.read_text()))
    # The scan body in ice_dfsn is a JAX implementation of a Fortran DO loop.
    functions.pop(('ice_dfsn', 'step'), None)
    # JAX cond bodies preserve Fortran fatal RETURN boundaries. Require their
    # direct use by lax.cond; this is not a blanket nested-helper exemption.
    for key in (('advance_microphys', 'advance_cloud_number'),
                ('advance_microphys', 'accumulate_statistics'),
                ('advance_ncm', 'finish_cloud_number'),
                ('microphys_solve', 'accumulate_implicit_statistics')):
        if key in functions:
            parent = functions[key[:-1]]
            assert any(isinstance(n, ast.Call) and ast.unparse(n.func) == 'jax.lax.cond'
                       and any(isinstance(a, ast.Name) and a.id == key[-1] for a in n.args)
                       for n in ast.walk(parent))
            functions.pop(key)
    assert list(functions) == list(routines)
    for key, function in functions.items():
        args = [arg.arg.lower() for arg in function.args.args]
        expected, intents = routines[key]
        assert args == [arg for arg in expected if arg in args], key
        assert not [arg for arg in expected if intents.get(arg) in ('in','inout') and arg not in args], key
        assert function.args.vararg is None and function.args.kwarg is None
        # Current source has no optional input/inout arguments. Its sole
        # optional PDF output is returned instead of becoming a Python input.
        assert not function.args.defaults and not function.args.kwonlyargs, key


@pytest.mark.parametrize('module,routine,expected', [
    ('morrison_microphys_module','morrison_microphys_driver',
     ['stats','hydromet_mc','hydromet_vel_zt','Ncm_mc','rcm_mc','rvm_mc','thlm_mc','rrm_auto_diag','rrm_accr_diag','rrm_evap_diag','Nrm_auto_diag','Nrm_evap_diag']),
    ('microphys_driver','calc_microphys_scheme_tendcies',
     ['stats','Nccnm','hydromet_mc','Ncm_mc','rcm_mc','rvm_mc','thlm_mc','hydromet_vel_zt','hydromet_vel_covar_zt_impc','hydromet_vel_covar_zt_expc','wprtp_mc','wpthlp_mc','rtp2_mc','thlp2_mc','rtpthlp_mc','Skw_zm_smooth']),
    ('microphys_init_cleanup','init_microphys',
     ['hydromet_dim','pdf_dim','hm_metadata','silhs_config_flags','vert_decorr_coef_out','corr_array_n_cloud','corr_array_n_below']),
])
def test_returned_state_order(module,routine,expected):
    tree=ast.parse((ROOT / f'clubb_jax/src/Microphys/{module}.py').read_text())
    function=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==routine)
    returned=[n for n in function.body if isinstance(n,ast.Return)][-1]
    assert [n.id for n in returned.value.elts] == expected


def test_driver_critical_order_and_no_obsolete_steps():
    code=(ROOT/'clubb_jax/src/advance_clubb_to_end.py').read_text()
    assert code.index('_advance_microphysics(state,') < code.index('_advance_radiation(state=',code.index('_advance_microphysics(state,'))
    section=code[code.index('def _advance_microphysics'):]
    assert section.index('_pdf_hydromet_microphys_prep_jit(') < section.index('_calc_microphys_scheme_tendcies_jit(') < section.index('_advance_microphys_jit(')
    assert 'kk_microphys_step' not in code and 'morrison_microphys_step' not in code


def test_runtime_interface_has_no_fortran_fallback():
    for source in MODULES:
        code=(ROOT/'clubb_jax'/source.relative_to(ROOT).with_suffix('.py')).read_text()
        assert 'import clubb_api' not in code
        assert 'pure_callback' not in code


@pytest.mark.parametrize('module', ['KK_microphys_module', 'morrison_microphys_module'])
def test_scheme_statistics_inventory_and_order(module):
    source=(ROOT/f'src/Microphys/{module}.F90').read_text()
    target=(ROOT/f'clubb_jax/src/Microphys/{module}.py').read_text()
    names=r'''\s*\(\s*["'](.*?)["']'''
    expected=re.findall(r'(?:stats_update|update_microphys_stat(?:_sfc)?)'+names,source,re.I)
    actual=re.findall(r'(?:stats\.update|update_microphys_stat(?:_sfc)?)'+names,target)
    assert actual == expected


@pytest.mark.parametrize('module,routine', [
    ('fill_holes', 'fill_holes_driver_api'),
    ('stats_clubb_utilities', 'stats_accumulate_hydromet_api'),
])
def test_updated_shared_hydrometeor_interfaces(module,routine):
    expected,intents=fortran_routines(ROOT/f'src/CLUBB_core/{module}.F90')[(routine,)]
    tree=ast.parse((ROOT/f'clubb_jax/src/CLUBB_core/{module}.py').read_text())
    function=python_routines(tree)[(routine,)]
    assert [a.arg.lower() for a in function.args.args] == [a for a in expected if intents.get(a) != 'out']


def test_fatal_dump_uses_actual_model_time():
    tree = ast.parse((ROOT/'clubb_jax/src/advance_clubb_to_end.py').read_text())
    call = next(n for n in ast.walk(tree) if isinstance(n, ast.Call)
                and isinstance(n.func, ast.Name) and n.func.id == 'write_adv_micro_errors')
    # The kernel receives a static startup-state time to avoid recompilation;
    # source diagnostics must report the actual failing model time instead.
    assert isinstance(call.args[3], ast.Name) and call.args[3].id == 'time_current'
