from scipy.constants import g
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

D = .075
W = .15
L = .5842
# L =.4342
# m = 1.036
m = 1
ds = .001
cg = np.array([0, 0, .0135])
# 120: cg = np.array([0, 0, .0242])
# 140: cg = np.array([0, 0, 00000])

X, Y, Z = np.meshgrid(np.arange(-L / 2, L / 2, ds),
                      np.arange(-W / 2, W / 2, ds),
                      np.arange(-.04, D, ds))


def hull(x, y):
    return .012 * np.tan(np.abs(39 * y) - 1.5) + .01


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
    waterMass = np.sum(waterMasses)
    return np.array([np.sum(X * waterMasses) / waterMass,
                     np.sum(Y * waterMasses) / waterMass,
                     np.sum(Z * waterMasses) / waterMass])


def rightingMoment(theta):
    '''Returns the righting moment'''
    rotation = np.reshape((np.array([1, 0, 0, 0, np.cos(
        theta), -np.sin(theta), 0, np.sin(theta), np.cos(theta)])), (3, 3))
    cbs = cb(theta)
    r = cbs - cg
    fb = np.dot(rotation, np.array([0, 0, m * g]))
    return (np.cross(r, fb)[0])


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

# print(rightingMoment(130 * np.pi / 180))

# newX, newY = np.meshgrid(np.arange(-L / 2, L / 2, ds),
#                       np.arange(-W / 2, W / 2, ds))

# newNewX = newX * hullMat[:,:,-1]
# newNewY = newY * hullMat[:,:,-1]

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(newNewX, newNewY, hull(newNewX, newNewY), cmap = cm.coolwarm)

# ax.figure.savefig('curve.png')

# plt.show()
waterMat = tilt(0)
plt.imshow(waterMat[:,0,:])
plt.show()
