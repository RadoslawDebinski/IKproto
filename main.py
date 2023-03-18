import numpy as np


# D-H parameters
# +------+---------+--------------------+
# |  d   |    a    |       alpha        |
# +------+---------+--------------------+
# |  0   |    0    |      pi / 2        |
# |  0   |    L1   |         0          |
# |  0   |    L2   |         0          |
# | -L3  |    L4   |    - pi / 2        |
# +------+---------+--------------------+

DHPARAMS = np.array([[0,       0,    np.pi / 2],
                     [0,    1450,            0],
                     [0,    1200,            0],
                     [-150,  310,   -np.pi / 2]])

def rotX(angle):
    rotMat = np.array([[1,             0,              0],
                       [0, np.cos(angle), -np.sin(angle)],
                       [0, np.sin(angle),  np.cos(angle)]])
    return rotMat


def rotY(angle):
    rotMat = np.array([[np.cos(angle),  0, np.sin(angle)],
                       [0,              1,             0],
                       [-np.sin(angle), 0, np.cos(angle)]])
    return rotMat


def rotZ(angle):
    rotMat = np.array([[np.cos(angle), -np.sin(angle), 0],
                       [np.sin(angle),  np.cos(angle), 0],
                       [0,                          0, 1]])
    return rotMat


# Forward coordinate system conversion - fCsc
def fCsc(x, y, z, alpha, beta, gamma):
    rot = rotX(alpha).dot(rotY(beta)).dot(rotZ(gamma))  # rotZ(gamma), rotY(beta)), rotX(alpha)
    rot = np.append(np.append(rot, [[0], [0], [0]], axis=1), [[0, 0, 0, 1]], axis=0)
    tran = np.array([[1, 0, 0, x],
                     [0, 1, 0, y],
                     [0, 0, 1, z],
                     [0, 0, 0, 1]])
    end = np.matmul(tran, rot)

    return end


# Inverse coordinate system conversion - iCsc
def iCsc(x, y, z, alpha, beta, gamma):
    rot = rotZ(-gamma).dot(rotY(-beta)).dot(rotX(-alpha))
    rot = np.append(np.append(rot, [[0], [0], [0]], axis=1), [[0, 0, 0, 1]], axis=0)
    tran = np.array([[1, 0, 0, -x],
                     [0, 1, 0, -y],
                     [0, 0, 1, -z],
                     [0, 0, 0,  1]])
    end = np.matmul(rot, tran)

    return end


def calcValH(u, d, l, a):
    H = [[np.cos(u),    -np.sin(u) * np.cos(a),      np.sin(u) * np.sin(a),     l * np.cos(u)],
         [np.sin(u),     np.cos(u) * np.cos(a),     -np.cos(u) * np.sin(a),     l * np.sin(u)],
         [0,                         np.sin(a),                  np.cos(a),                 d],
         [0,                                 0,                          0,                 1]]

    return np.array(H)



if __name__ == '__main__':

    # INPUT DATA
    xb, yb, zb, rXb, rYb, rZb = 0, 0, 450, 0, 0, np.pi
    x4, y4, z4, rX4, rY4, rZ4 = -804.847, 303.021, 2438.528, np.pi, 0, np.pi/6

    h0toB = iCsc(xb, yb, zb, rXb, rYb, rZb)
    hBto4 = fCsc(x4, y4, z4, rX4, rY4, rZ4)

    q5 = rZ4

    h0to4 = h0toB.dot(hBto4)

    # pw DATA
    werY = np.array([[h0to4[0][1]], [h0to4[1][1]], [h0to4[2][1]]])
    werZ = np.array([[h0to4[0][2]], [h0to4[1][2]], [h0to4[2][2]]])
    tran = np.array([[h0to4[0][3]], [h0to4[1][3]], [h0to4[2][3]]])

    # pw COMPUTE and DH PARAMS
    L1 = DHPARAMS[1][1]
    L2 = DHPARAMS[2][1]
    L3 = DHPARAMS[3][0]
    L4 = DHPARAMS[3][1]

    tranWriY = L3*werY
    tranWriZ = L4*werZ

    pw = tran + tranWriY - tranWriZ

    # compute first 3 angles analytically
    xp, yp, zp = pw[0], pw[1], pw[2]
    r = np.sqrt(np.square(xp) + np.square(yp))
    s = np.sqrt(np.square(r) + np.square(zp))

    alpha = np.arctan2(zp, r)
    beta = np.arccos((np.square(s) + np.square(L1) - np.square(L2))/(2 * s * L1))
    gamma = np.arccos((np.square(L1) + np.square(L2) - np.square(s))/(2 * L1 * L2))

    q1 = (np.arctan2(yp, xp))[0]
    q2 = (alpha + beta)[0]
    q3 = (gamma - np.pi)[0]

    # compute 4th angle by dh matrices

    d1, d2, d3, d4 = DHPARAMS[:, 0]  # Delta
    l1, l2, l3, l4 = DHPARAMS[:, 1]  # Lambda
    a1, a2, a3, a4 = DHPARAMS[:, 2]  # Alpha

    h0to1 = calcValH(q1, d1, l1, a1)
    h1to2 = calcValH(q2, d2, l2, a2)
    h2to3 = calcValH(q3, d3, l3, a3)

    h0to3 = h0to1.dot(h1to2).dot(h2to3)

    h3to4 = np.linalg.inv(h0to3).dot(h0to4)

    sinq4 = h3to4[1][3]
    cosq4 = h3to4[0][3]

    tgq4 = sinq4 / cosq4

    q4 = np.arctan2(sinq4, cosq4)

    print(f"angle q1 = {np.round(q1 * 180 / np.pi, 5)}°")
    print(f"angle q2 = {np.round(q2 * 180 / np.pi, 5)}°")
    print(f"angle q3 = {np.round(q3 * 180 / np.pi, 5)}°")
    print(f"angle q4 = {np.round(q4 * 180 / np.pi, 5)}°")
    print(f"angle q5 = {np.round(q5 * 180 / np.pi, 5)}°")