import numpy as np


def normalize(p: np.ndarray) -> np.ndarray:
    return (p % 1.0) % 1.0


def equivalent_sets(setA, setB, **kwargs):
    pass


def list_to_str(l, sep='\t\t', decimals=5, width=8):
    out_str = ""
    for (index, val) in enumerate(l):
        if isinstance(val, int):
            out_str += f"{val:{width}}"
        else:
            out_str += f"{val:{width}.{decimals}f}"
        
        if index + 1 < len(l):
            out_str += sep
     
    return out_str


# Positions
def compute_unit_cell_position(position: np.ndarray) -> np.ndarray:
    return normalize(position)


def compute_shift_position(position: np.ndarray):
    unit_cell_position = compute_unit_cell_position(position=position)
    shift = np.round(position - unit_cell_position).astype(int)

    return (unit_cell_position, shift)


# Distances
def periodic_diff(p1: np.ndarray, p2: np.ndarray) -> np.ndarray:
    diff1 = normalize(p1 - p2)
    diff2 = normalize(p2 - p1)
    min_diff = np.minimum(diff1, diff2)

    return min_diff
