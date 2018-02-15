print(rightingMoment(125.50 * np.pi / 180))
from scipy.constants import g
import matplotlib.pyplot as plt
import numpy as np

# D = .05  # depth (Z)
# W = .075  # width (Y)
# L = .5842  # length (X)
# m = 1.418  # mass
# ds = .005
# cg = np.array([0, 0, 0.02])

# mvp
# D = .1
# W = .524
# L = .
# m = 1.4
# ds = .003
# cg = np.array([0, 0, .055])

# D = .05
# W = .1
# L = .15
# m = .7
# ds = .003
# cg = np.array([0,0,.025])

# real
# D = .10
# W = .1524
# L = .58
# m = .8997
# ds = .003
# cg = np.array([0, 0, .035])



# narrow real
D = .06
W = .11
L = .58
m = 1.418
ds = .003
cg = np.array([0, 0, .025])

X, Y, Z = np.meshgrid(np.arange(-L / 2, L / 2, ds),
                      np.arange(-W / 2, W / 2, ds),
                      np.arange(0, D, ds))


def hull(x, y):
    # return -(W**2 / 4 - y ** 2) ** .5 + W / 2
    # return D * (4 * x ** 2 / L ** 2 + 4 * y ** 2 / W ** 2)
    return - ((W / 2 - 128 * W * x ** 8 / L ** 8) ** 2 - y ** 2) ** .5 + W / 2


hullMat = Z > hull(X, Y)

volume = np.sum(hullMat * ds ** 3)


def tilt(theta):
    b = intercept(theta)
    water = waterLine(theta, b, Y)
    if np.abs(theta) < np.pi / 2:
        underWater = Z < water
    else:
        underWater = Z > water
    return np.logical_and(underWater, hullMat)


def intercept(theta):
    def massDiff(b):
        if np.abs(theta) < np.pi / 2:
            underWater = Z < waterLine(theta, b, Y)
        else:
            underWater = Z > waterLine(theta, b, Y)
        hullWater = np.logical_and(underWater, hullMat)
        massWater = np.sum(hullWater) * ds ** 3 * 1000
        return np.abs(massWater - m)
    maxB = np.abs(np.tan(theta)) * W / 2 + D
    intercept = np.argmin([massDiff(b)
                           for b in np.arange(0, maxB, ds)]) * ds
    return intercept


def waterLine(theta, b, y):
    return np.tan(theta) * y + b


def cb(theta):
    '''Returns the center of bouyancy'''
    waterMat = tilt(theta)
    waterMasses = waterMat * ds ** 3 * 1000
    # plt.imshow(waterMat[18,:,:])
    # plt.imshow(hullMat[18,:,:])
    # plt.show()
    waterMass = np.sum(waterMasses)
    return np.array([np.sum(X * waterMasses) / waterMass,
                     np.sum(Y * waterMasses) / waterMass,
                     np.sum(Z * waterMasses) / waterMass])


def rightingMoment(theta):
    '''Returns the righting moment'''
    # cb = cb(theta)[1:]
    # r = (cb - cm[1:])
    # fb = np.array([0, 0, g * m])  # bouyant force
    rotation = np.reshape((np.array([1, 0, 0, 0, np.cos(
        theta), -np.sin(theta), 0, np.sin(theta), np.cos(theta)])), (3, 3))
    cbs = cb(theta)
    r = cbs - cg
    # print(cb)
    fb = np.dot(rotation, np.array([0, 0, m * g]))
    return (np.cross(r, fb)[0])
    # return [0,0,0]


# cb(0)

# def intercept(
# plt.imshow(tilt(np.pi / 3)[:, 100, :])
# plt.show()
# plt.imshow(hullMat[:, :, :])
# ax = plt.gca()
# ax.set_aspect('equal')
# plt.imshow(hullMat[:, :, 36])
# plt.show()
# print(np.sum((hullMat)))

angles = np.arange(0, 360, 4) * np.pi / 180
moments = [rightingMoment(angle) for angle in angles]
plt.plot(np.degrees(angles), moments)
plt.show()

# print(rightingMoment(125.368 * np.pi / 180))

# print(rightingMoment(131 * np.pi / 180))

# print(rightingMoment(5 * np.pi / 9))
