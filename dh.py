import numpy as np
import sympy as sym


# D-H parameters
# +------+---------+--------------------+
# |  d   |    l    |       alpha        |
# +------+---------+--------------------+
# |  0   |    0    |      pi / 2        |
# |  0   |    L1   |         0          |
# |  0   |    L2   |         0          |
# | -L3  |    L4   |    - pi / 2        |
# +------+---------+--------------------+


def calcSymH(u, d, l, a):
    H = [[sym.cos(u),   -sym.sin(u) * np.cos(a),    sym.sin(u) * np.sin(a),     l * sym.cos(u)],
         [sym.sin(u),   sym.cos(u) * np.cos(a),     -sym.cos(u) * np.sin(a),    l * sym.sin(u)],
         [0,            np.sin(a),                  np.cos(a),                  d],
         [0,            0,                          0,                          1]]

    return np.array(H)


def calcValH(u, d, l, a):
    H = [[np.cos(u),   -np.sin(u) * np.cos(a),     np.sin(u) * np.sin(a),      l * np.cos(u)],
         [np.sin(u),   np.cos(u) * np.cos(a),      -np.cos(u) * np.sin(a),     l * np.sin(u)],
         [0,            np.sin(a),                   np.cos(a),                   d],
         [0,            0,                            0,                            1]]

    return np.array(H)


if __name__ == '__main__':
    # Robot Lengths
    L1 = 1450
    L2 = 1200
    L3 = 150
    L4 = 200

    # Upsilon
    u1 = sym.Symbol('u1')
    u2 = sym.Symbol('u2')
    u3 = sym.Symbol('u3')
    u4 = sym.Symbol('u4')

    # Delta
    d1 = 0
    d2 = 0
    d3 = 0
    d4 = -L3

    # Lambda
    l1 = 0
    l2 = L1
    l3 = L2
    l4 = L4

    # Alpha
    a1 = np.pi / 2
    a2 = 0
    a3 = 0
    a4 = -np.pi / 2

    # -> D-H Transformation Matrix - Symbolic Calculations
    H1 = calcSymH(u1, d1, l1, a1)
    H2 = calcSymH(u2, d2, l2, a2)
    H3 = calcSymH(u3, d3, l3, a3)
    H4 = calcSymH(u4, d4, l4, a4)

    H_sym = sym.simplify(H1.dot(H2).dot(H3).dot(H4))
    print(f"Symbolic H: \n{H_sym}")

    # -> D-H Transformation Matrix - Checking Values
    H1 = calcValH(np.pi/2, d1, l1, a1)
    H2 = calcValH(np.pi, d2, l2, a2)
    H3 = calcValH(0, d3, l3, a3)
    H4 = calcValH(np.pi, d4, l4, a4)

    H_val = H1.dot(H2).dot(H3).dot(H4)
    detH = np.linalg.det(H_val)
    # print(f"H with values: \n{H_val}")
    # print(f"det(H) = {detH}")
