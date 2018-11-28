import numpy as np
import cmath
import re

def normalize(x):
    return x/np.linalg.norm(x)


def get_trailing_numbers(s, zero=False):
    m = re.search(r'\d+$', s)
    if m:
        return int(m.group())
    elif zero:
        return 0
    else:
        return None

