"""Read normal-space correlation matrices (corr_varnce_input_reader.F90).

File parsing is host-only; indices are zero based.
"""
from pathlib import Path
import numpy as np
from clubb_jax.src.CLUBB_core.corr_varnce_module import get_corr_var_index


def read_correlation_matrix(iunit, input_file, pdf_dim, hm_metadata, corr_array_n):
    lines = [line.split('!')[0].strip() for line in Path(input_file).read_text().splitlines()]
    lines = [line for line in lines if line]
    names = lines[0].split()
    values = np.array([[float(x) for x in line.split()] for line in lines[1:]])
    if values.shape != (len(names), len(names)):
        raise ValueError(f"Correlation matrix must have equal rows and columns: {input_file}")
    corr_array_n = np.eye(pdf_dim)
    for i in range(1, len(names)):
        var_index1 = get_corr_var_index(names[i], hm_metadata)
        if var_index1 >= 0:
            for j in range(i):
                var_index2 = get_corr_var_index(names[j], hm_metadata)
                if var_index2 >= 0:
                    corr_array_n[max(var_index1, var_index2), min(var_index1, var_index2)] = values[j, i]
    return corr_array_n
