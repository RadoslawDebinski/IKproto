import numpy as np
import sympy as sym
from prettytable import PrettyTable

# D-H parameters
# +------+---------+--------------------+
# |  d   |    l    |       alpha        |
# +------+---------+--------------------+
# |  0   |    0    |      pi / 2        |
# |  0   |    L1   |         0          |
# |  0   |    L2   |         0          |
# | -L3  |    L4   |    - pi / 2        |
# +------+---------+--------------------+

DHPARAMS = np.array([[0,       0,   np.pi / 2],
                     [0,    1450,           0],
                     [0,    1200,           0],
                     [-150,  310,  -np.pi / 2]])


def calcSymH(u, d, l, a):
    H = [[sym.cos(u),   -sym.sin(u) * np.cos(a),     sym.sin(u) * np.sin(a),    l * sym.cos(u)],
         [sym.sin(u),    sym.cos(u) * np.cos(a),    -sym.cos(u) * np.sin(a),    l * sym.sin(u)],
         [0,                          np.sin(a),                  np.cos(a),                 d],
         [0,                                  0,                          0,                 1]]

    return np.array(H)


def calcValH(u, d, l, a):
    H = [[np.cos(u),    -np.sin(u) * np.cos(a),      np.sin(u) * np.sin(a),     l * np.cos(u)],
         [np.sin(u),     np.cos(u) * np.cos(a),     -np.cos(u) * np.sin(a),     l * np.sin(u)],
         [0,                         np.sin(a),                  np.cos(a),                 d],
         [0,                                 0,                          0,                 1]]

    return np.array(H)


if __name__ == '__main__':
    # Upsilon
    u1 = sym.Symbol('u1')
    u2 = sym.Symbol('u2')
    u3 = sym.Symbol('u3')
    u4 = sym.Symbol('u4')

    d1, d2, d3, d4 = DHPARAMS[:, 0]  # Delta
    l1, l2, l3, l4 = DHPARAMS[:, 1]  # Lambda
    a1, a2, a3, a4 = DHPARAMS[:, 2]  # Alpha

    # -> D-H Transformation Matrix - Symbolic Calculations
    H1 = calcSymH(u1, d1, l1, a1)
    H2 = calcSymH(u2, d2, l2, a2)
    H3 = calcSymH(u3, d3, l3, a3)
    H4 = calcSymH(u4, d4, l4, a4)

    H_sym = sym.simplify(H1.dot(H2).dot(H3).dot(H4))
    H_tab = PrettyTable()
    H_tab.add_column("1st", H_sym[:, 0])
    H_tab.add_column("2nd", H_sym[:, 1])
    H_tab.add_column("3rd", H_sym[:, 2])
    H_tab.add_column("4th", H_sym[:, 3])
    print("Symbolic H: \n", H_tab)

    # -> D-H Transformation Matrix - Checking Values
    H1 = calcValH(np.pi / 2, d1, l1, a1)
    H2 = calcValH(np.pi, d2, l2, a2)
    H3 = calcValH(0, d3, l3, a3)
    H4 = calcValH(np.pi, d4, l4, a4)
    H_val = H1.dot(H2).dot(H3).dot(H4)

    # detH = np.linalg.det(H_val)
    # print(f"H with values: \n{H_val}")
    # print(f"det(H) = {detH}")
