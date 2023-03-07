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


if __name__ == '__main__':
    # h0to5 = np.matmul(iCsc(0, 0, 450, 0, 0, np.pi), fCsc(-1925, -140, 1395.737, 0, np.pi, 0))

    h0toB = iCsc(0, 0, 450, 0, 0, np.pi)
    hBto5 = fCsc(-1925, -140, 1395.737, 0, np.pi, 0)

    h0to5 = hBto5 # h0toB.dot(hBto5)

    # # # TESTS STEP 1
    # # # SIMPLE TEST FOR "Forward coordinate system conversion - fCsc" - PASS!
    # # print(np.matmul(fCsc(0, 0, 10, 0, np.pi, np.pi), [[1], [1], [1], [1]]))
    # print(np.matmul(fCsc(2, 0, 0, 0, np.pi/2, np.pi/2), [[1], [2], [3], [1]]))
    #
    # # # TESTS STEP 2
    # # # SIMPLE TEST FOR "Inverse coordinate system conversion - iCsc" - PASS!
    # # print(np.matmul(iCsc(0, 0, 10, 0, np.pi, np.pi), [[1], [1], [1], [1]]))
    # print(np.matmul(iCsc(2, 0, 0, 0, np.pi/2, np.pi/2), [[1], [2], [3], [1]]))

    # INPUT DATA
    # h0toB(0, 0, 450, 0, 0, np.pi)
    # hBto5(-1925, -140, 1395.737, 0, np.pi, 0)
    # OUTPUT DATA
    # -1925, 10, 1705.737

    # pw DATA
    werX = np.array([[h0to5[0][0]], [h0to5[1][0]], [h0to5[2][0]]])
    werZ = np.array([[h0to5[0][2]], [h0to5[1][2]], [h0to5[2][2]]])
    tran = np.array([[h0to5[0][3]], [h0to5[1][3]], [h0to5[2][3]]])
    # # OUTPUT pw DATA
    # # print(tran)
    # # print(werX)
    # # print(werZ)
    # # print(DHPARAMS)
    #
    # # pw COMPUTE and DH PARAMS
    L3 = -DHPARAMS[3][0]
    L4 = DHPARAMS[3][1]
    print(f"{L3} {L4}")
    pw = tran - (L4*werX) - ((-L3)*werZ)
    print(pw)