import numpy as np

mass = 1.43
length = .5841
# density = 720
width = 10.16
ds = .001
height = 12.16

X, Y, Z = np.meshgrid(np.arange(- length / 2, length / 2, ds),
                      np.arange(- width / 2, width / 2, ds), np.arange(0, height, ds))

hull = (- (width ** 2 / 4 - Y ** 2) ** .5 + width / 2) < Z

volume = np.sum(hull) * ds ** 3

density = mass / volume

masses = hull * ds ** 3 * density

cg = np.array([0, 0, np.sum(masses * Z) / mass])

theta = 2 * np.pi / 9
