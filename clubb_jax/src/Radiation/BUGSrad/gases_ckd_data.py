"""Parse the BUGSrad correlated-k coefficient tables from gases_ckd_data.h.

The Fortran data file (~1522 lines) is too large to hand-transcribe reliably, so this module PARSES
the Fortran `data` statements at import and builds the coefficient arrays. Three statement forms:
  - 1D:  `data var / v1, v2, ... /`                                   (hk1..hk18, fk1o3)
  - 3D:  `data((var(T_k,j,i),i=1,KG),j=1,NUMP) / ... /`  (per T-coef)  (c*h2o, c*o3 — i fastest)
  - 2D:  `data((var(k,j),j=1,NUMP),k=1,NUMT) / ... /`                  (c10ch4, c10n2o — j fastest)
Each 3D array `var(NUMT=3, NUMP, KG)` is filled by 3 statements (T_1/T_2/T_3 = the a/b/c poly coefs).
Validated by shape + value-count + spot checks in tests/test_bugsrad.py.
"""
import os
import re
import numpy as np

# parameters (gases_ckd_data.h:16-53)
MB, MBS, MBIR = 18, 6, 12
KG = [10, 8, 12, 7, 12, 5, 2, 3, 4, 4, 3, 5, 2, 10, 12, 7, 7, 8]   # k-samples per band
NUMPS, NUMTS, NUMPIR, NUMTIR = 11, 3, 19, 3
STANPS = np.array([10.0, 15.8, 25.1, 39.8, 63.1, 100.0, 158.0, 251.0, 398.0, 631.0, 1000.0])
STANPIR = np.array([0.251, 0.398, 0.631, 1.000, 1.58, 2.51, 3.98, 6.31, 10.0, 15.8, 25.1, 39.8,
                    63.1, 100.0, 158.0, 251.0, 398.0, 631.0, 1000.0])

# symbol table for the implied-DO bounds
_SYM = {'NUMPS': NUMPS, 'NUMTS': NUMTS, 'NUMPIR': NUMPIR, 'NUMTIR': NUMTIR,
        'T_1': 1, 'T_2': 2, 'T_3': 3, **{f'KG_{k + 1}': KG[k] for k in range(MB)}}

_DATA_H = os.path.join(os.path.dirname(__file__), '..', '..', '..', '..',
                       'clubb_release', 'src', 'Radiation', 'BUGSrad', 'gases_ckd_data.h')


def _resolve(tok):
    tok = tok.strip()
    return int(tok) if re.fullmatch(r'\d+', tok) else _SYM[tok]


def _logical_lines(path):
    """Strip comments + #directives, join `&`-continued lines into logical statements."""
    out, buf = [], ''
    for raw in open(path):
        line = raw.split('!', 1)[0].rstrip('\n')
        if line.lstrip().startswith('#'):
            continue
        line = line.rstrip()
        if line.endswith('&'):
            buf += ' ' + line[:-1]
        else:
            buf += ' ' + line
            if buf.strip():
                out.append(buf.strip())
            buf = ''
    if buf.strip():
        out.append(buf.strip())
    return out


def _parse_values(s):
    return np.array([float(v) for v in s.split(',') if v.strip()], dtype=np.float64)


def load_tables():
    """Return a dict {name: ndarray} of the correlated-k tables (1D hk*/fk1o3, plus the c* arrays:
    3D c*h2o/c*o3 of shape (3, NUMP, KG_band), 2D c*ch4/c*n2o of shape (NUMTIR, NUMPIR))."""
    tabs = {}
    for stmt in _logical_lines(os.path.normpath(_DATA_H)):
        m = re.match(r'data\s*(.*)$', stmt)
        if not m or not (m.group(1).startswith('((') or re.match(r'[a-z]', m.group(1))):
            continue
        body = m.group(1)
        # split header / values / (trailing)
        parts = body.split('/')
        header = parts[0].strip()
        values = _parse_values(parts[1]) if len(parts) > 1 else None
        if values is None:
            continue
        if header.startswith('(('):
            # implied-DO. Extract var + the (Tk,j,i) or (k,j) + the two loop specs.
            mm = re.match(r'\(\(\s*([a-z0-9_]+)\s*\(([^)]*)\)\s*,\s*([A-Za-z]+)\s*=\s*\d+\s*,'
                          r'\s*([A-Za-z0-9_]+)\s*\)\s*,\s*([A-Za-z]+)\s*=\s*\d+\s*,\s*([A-Za-z0-9_]+)\s*\)',
                          header)
            var, idxspec, inner_v, inner_hi, outer_v, outer_hi = mm.groups()
            inner_n, outer_n = _resolve(inner_hi), _resolve(outer_hi)
            ids = [t.strip() for t in idxspec.split(',')]
            if ids[0].startswith('T_'):                       # 3D: var(T_k, j, i), i inner, j outer
                t = _SYM[ids[0]] - 1
                arr = tabs.setdefault(var, np.zeros((NUMTS if outer_n == NUMPS else NUMTIR,
                                                     outer_n, inner_n)))
                arr[t] = values.reshape(outer_n, inner_n)     # [j (outer), i (inner)]
            else:                                             # 2D: var(k, j), j inner, k outer
                tabs[var] = values.reshape(outer_n, inner_n)  # [k (outer), j (inner)]
        else:
            tabs[header.strip()] = values                     # 1D
    return tabs


_TABLES = None


def tables():
    """Cached parse of gases_ckd_data.h."""
    global _TABLES
    if _TABLES is None:
        _TABLES = load_tables()
    return _TABLES
