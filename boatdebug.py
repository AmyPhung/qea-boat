from scipy.constants import g
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# mvp
# D = .1
# W = .524
# L = .
# m = 1.4
# ds = .003
# cg = np.array([0, 0, .055])

# real
# D = .10
# W = .1524
# L = .58
# m = .8997
# ds = .003
# cg = np.array([0, 0, .035])

# narrow real
# D = .075
# W = .15
# L = .7112
# m = 1
# ds = .002
# cg = np.array([0, 0, .0170])

# D = 1
# W = 1
# L = 1
# m = 500
# ds = 0.015
# cg = np.array([0,0,.5])

D = .075
W = .15
L = .5842
m = 1
ds = .003
# 130: cg = np.array([0, 0, .01070])
# cg = np.array([0, 0, .023])
# cg = np.array([0, 0, -.001])

X, Y, Z = np.meshgrid(np.arange(-L / 2, L / 2, ds),
                      np.arange(-W / 2, W / 2, ds),
                      np.arange(-.04, D, ds))


def hull(x, y):
    # return -(W**2 / 4 - y ** 2) ** .5 + W / 2
    # return D * (4 * x ** 2 / L ** 2 + 4 * y ** 2 / W ** 2)
    # return - ((W / 2 - 128 * W * x ** 8 / L ** 8) ** 2 - y ** 2) ** .5 + W / 2
    # r = W / 2
    # r = 0.075
    # return 0.05 * np.tan(np.abs(22 * y) - .9) + .01
    return .018 * np.tan(np.abs(33 * y) - 1.2) + .005
    # return y
    # return 75 * np.tan(np.abs(y) / .075 * 1.2 - .6) + .01
    # return np.tan(y)
    # return y

    # cube
    # return 0


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
# plt.imshow(hullMat[:, 58, :])
# ax = plt.gca()
# ax.set_aspect('equal')
# plt.imshow(tilt(0)[:, 356, :])
# plt.show()
# print(np.sum((hullMat)))

# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')
# angles = np.arange(0, 180, 1.1) * np.pi / 180
# ax = plt.gca()
# moments = [rightingMoment(angle) for angle in angles]
# plt.plot(np.degrees(angles), moments)
# ax.set_xlabel(r'Angle (\circ)')
# ax.set_ylabel(r'Righting Moment ($N \cdot m$)')
# plt.show()
# ax.figure.savefig('curve.png')

print(rightingMoment(140 * np.pi / 180))

# newX, newY = np.meshgrid(np.arange(-L / 2, L / 2, ds),
#                       np.arange(-W / 2, W / 2, ds))

# newNewX = newX * hullMat[:,:,-1]
# newNewY = newY * hullMat[:,:,-1]

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(newNewX, newNewY, hull(newNewX, newNewY), cmap = cm.coolwarm)

# ax.figure.savefig('badlyscaledbutbiggerhull.png')

# plt.show()
