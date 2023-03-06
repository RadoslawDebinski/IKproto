import numpy as np


# D-H parameters
# +------+---------+--------------------+
# |  d   |    a    |       alpha        |
# +------+---------+--------------------+
# |  0   |    0    |      pi / 2        |
# |  0   |    L1   |         0          |
# |  0   |    L1   |         0          |
# | -L3  |    L4   |    - pi / 2        |
# +------+---------+--------------------+

DHPARAMS = np.array([[0,       0,    np.pi / 2],
                     [0,    1450,            0],
                     [0,    1200,            0],
                     [-150,  200,   -np.pi / 2]])

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


def hBto5(x, y, z, alpha5, beta5, gamma5):
    rot = np.matmul(np.matmul(rotX(alpha5), rotY(beta5)), rotZ(gamma5))
    rotAndTran = np.append(np.append(rot, [[x], [y], [z]], axis=1), [[0, 0, 0, 1]], axis=0)
    return rotAndTran


def h0toB(x, y, z, alpha0, beta0, gamma0):
    rotAndTran = hBto5(-x, -y, -z, -alpha0, -beta0, -gamma0)
    return rotAndTran


if __name__ == '__main__':
    h0to5 = np.matmul(h0toB(1, 1, 1, 0, 0, 0), hBto5(1, 2, 3, np.pi / 2, np.pi / 2, 0))
    print(h0to5)
    werX = np.array([[h0to5[0][0]], [h0to5[1][0]], [h0to5[2][0]]])
    werZ = np.array([[h0to5[0][2]], [h0to5[1][2]], [h0to5[2][2]]])
    tran = np.array([[h0to5[0][3]], [h0to5[1][3]], [h0to5[2][3]]])

    print(tran)
    print(werX)
    print(werZ)
    print(DHPARAMS)
    L3 = -DHPARAMS[3][0]
    L4 = DHPARAMS[3][1]
    print(f"{L3} {L4}")

    pw = tran - (L4*werX) - ((-L3)*werZ)

    print(pw)