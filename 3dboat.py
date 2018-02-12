import numpy as np
from scipy.constants import g


class boat:
    '''Object for a boat. All units are SI units, and angles are in radians, except avs.'''

    def __init__(self, D, L, W, m):
        self.D = D  # depth (Z)
        self.W = W  # width (Y)
        self.L = L  # length (X)
        self.m = m  # mass

        # Precision value
        self.ds = 0.001

        # meshgrids
        X, Y, Z = np.meshgrid(np.arange(-W / 2, W / 2, self.ds),
                              np.arange(-L / 2, L / 2, self.ds),
                              np.arange(0, D, self.ds))

        # Logical matrix of the hull
        self.hullMat = Z > self.hull(X, Y)

        masses = self.hullMat * m / self.volume

        self.volume = np.sum(self.hullMat) * self.ds ** 3

        # Position of the center of mass.
        self.cm = np.array([0,
                            0,
                            np.sum(Z * masses) / m])

    def hull(self, x, y):
        '''Returns the shape of the hull'''
        return self.b * ((2 * x / self.L) ** 2 + (2 * y / self.W) ** 2)

    def tilt(self, theta):
        '''Returns logical matrix of the part of the hull under
        water given a tilt'''
        b = self.intercept(theta)
        water = self.waterLine(theta, b, self.Y)
        if np.abs(theta) < np.pi / 2:
            underWater = self.Z < water
        else:
            underWater = self.Z > water
        return np.logical_and(underWater, self.hullMat)

    def cb(self, theta):
        '''Returns the center of bouyancy'''
        waterMat = self.tilt(theta)
        waterMasses = waterMat * self.ds ** 2 * self.L * 1000
        waterMass = 1000
        return np.array([0,
                         0,
                         np.sum(self.Z * waterMasses) / waterMass])

    def avs(self):
        '''returns avs in degrees. Is a method and not a field because this will
        probably take a while to calculate'''
        return np.argmin([rightingMoment(angle * np.pi / 180) for angle in np.arange(0, 90)])

    def rightingMoment(self, theta):
        '''Returns the righting moment'''
        r = self.cb(theta) - self.cm
        fb = np.array([0, 0, g * self.m])  # bouyant force
        return np.cross(r, fb)

    def intercept(self, theta):
        '''tilt the boat and return the y - intercept of the
        water'''

        # calculate the mass difference between hull and water
        def massDiff(b):
            if np.abs(theta) < np.pi / 2:
                underWater = self.Z < self.waterLine(theta, b, self.Y)
            else:
                underWater = self.Z > self.waterLine(theta, b, self.Y)
            hullWater = np.logical_and(underWater, self.hullMat)
            massWater = np.sum(hullWater) * self.ds ** 2 * self.L * 1000
            return np.abs(massWater - self.m)

        # maximum intercept to test, which is the one which
        # crosses one vertice of the boat.
        maxB = np.abs(np.tan(theta)) * self.W / 2 + self.D

        # Find the intercept with the smallest mass diff.
        intercept = np.argmin([massDiff(b)
                               for b in np.arange(0, maxB, self.ds)]) * self.ds
        return intercept

    # Waterline function (YZ plane)
    def waterLine(self, theta, b, y):
        '''Return the z-value of the waterline.'''
        return np.tan(theta) * y + b
