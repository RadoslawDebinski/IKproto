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


if __name__ == '__main__':
    # h0to5 = np.matmul(iCsc(0, 0, 450, 0, 0, np.pi), fCsc(-1925, -140, 1395.737, 0, np.pi, 0))

    h0toB = iCsc(0, 0, 450, 0, 0, np.pi)
    hBto5 = fCsc(-804.847, 303.021, 2438.528, np.pi, 0, np.pi/6)
    q5 = np.array([np.pi/6])
    h0to5 = h0toB.dot(hBto5)

    # # # TESTS STEP 1
    # # # SIMPLE TEST FOR "Forward coordinate system conversion - fCsc" - PASS!
    # print(np.matmul(fCsc(0, 0, 10, 0, np.pi, np.pi), [[1], [1], [1], [1]]))
    # print(np.matmul(fCsc(2, 0, 0, 0, np.pi/2, np.pi/2), [[1], [2], [3], [1]]))
    #
    # # # TESTS STEP 2
    # # # SIMPLE TEST FOR "Inverse coordinate system conversion - iCsc"
    # print(np.matmul(iCsc(0, 0, 10, 0, np.pi, np.pi), [[1], [1], [1], [1]]))
    # print(np.matmul(iCsc(2, 0, 0, 0, np.pi/2, np.pi/2), [[1], [2], [3], [1]]))

    # EXAMPLE OF INPUT DATA 1
    # h0toB(0, 0, 450, 0, 0, np.pi)
    # hBto5(-1925, -140, 1395.737, np.pi, 0, 0)
    # EXAMPLE OF OUTPUT DATA 1
    # pw = -1925, 10, 1705.737 - PASS!

    # EXAMPLE OF INPUT DATA 2
    # h0toB(0, 0, 0, 0, 0, 0)
    # hBto5(-2194.111, -140, 2016.32, np.pi, -np.pi / 2, 0)
    # EXAMPLE OF OUTPUT DATA 2
    # pw = -1884.111, 10, 2016.32 - PASS!

    # EXAMPLE OF INPUT DATA 3
    # h0toB(0, 0, 0, 0, 0, 0)
    # hBto5(-804.847, 303.021, 2438.528, np.pi, 0, np.pi/6)
    # EXAMPLE OF OUTPUT DATA 3
    # pw = -729.847, 432.924, 2748.528 - PASS!

    # EXAMPLE OF INPUT DATA 4
    # h0toB(0, 0, 0, 0, 0, 0)
    # hBto5(-140, 1764.23, 795.737, np.pi, 0, np.pi/2)
    # EXAMPLE OF OUTPUT DATA 4
    # pw = 10, 1764.23, 1105.737 - PASS!

    # EXAMPLE OF INPUT DATA 5
    # h0toB(0, 0, 0, 0, 0, 0)
    # hBto5(-1057.125, -1255.115, 666.506, -np.pi/2, -np.pi / 2 -np.pi / 4, np.pi/2)
    # EXAMPLE OF OUTPUT DATA 5
    # pw = -943.988, -929.845, 666.506 - PASS!

    # EXAMPLE OF INPUT DATA 6
    # h0toB(0, 0, 0, 0, 0, 0)
    # hBto5(-1985.483, -676.947, 2543.036, -np.pi - 1 * np.pi / 12, -np.pi / 2 - np.pi / 4, np.pi/2 + 7 * np.pi / 12)
    # EXAMPLE OF OUTPUT DATA 6
    # pw = -1812.572, -475.324, 2323.833

    # pw DATA
    werY = np.array([[h0to5[0][1]], [h0to5[1][1]], [h0to5[2][1]]])
    werZ = np.array([[h0to5[0][2]], [h0to5[1][2]], [h0to5[2][2]]])
    tran = np.array([[h0to5[0][3]], [h0to5[1][3]], [h0to5[2][3]]])

    # # pw COMPUTE and DH PARAMS
    L1 = DHPARAMS[1][1]
    L2 = DHPARAMS[2][1]
    L3 = DHPARAMS[3][0]
    L4 = DHPARAMS[3][1]

    tranWriY = L3*werY
    tranWriZ = L4*werZ

    pw = tran + tranWriY - tranWriZ

    xp, yp, zp = pw[0], pw[1], pw[2]
    r = np.sqrt(np.square(xp) + np.square(yp))
    s = np.sqrt(np.square(r) + np.square(zp))

    alpha = np.arctan2(zp, r)
    beta = np.arccos((np.square(s) + np.square(L1) - np.square(L2))/(2 * s * L1))
    gamma = np.arccos((np.square(L1) + np.square(L2) - np.square(s))/(2 * L1 * L2))

    q1 = np.arctan2(yp, xp)
    q2 = np.pi / 2 - alpha - beta
    q3 = np.pi - gamma

    print(f"angle q1 = {np.round(q1[0] * 180 / np.pi, 5)}°")
    print(f"angle q2 = {np.round(q2[0] * 180 / np.pi, 5)}°")
    print(f"angle q3 = {np.round(q3[0] * 180 / np.pi, 5)}°")
    # print(f"angle q4 = {np.round(180 / np.pi, 5)}°")
    print(f"angle q5 = {np.round(q5[0] * 180 / np.pi, 5)}°")



