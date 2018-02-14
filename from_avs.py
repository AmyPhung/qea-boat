import numpy as np
from scipy.constants import g
from scipy.integrate import quad

supportWidth = .01
length = .59  # approx 23.2 in
thickness = .003175  # .125 in
teakDensity = 720
ds = 0.001


def ellipse(w, x):
    # y-coordinate of ellips of width w and length length.
    return w / 2 / length * (length ** 2 - 4 * x ** 2) ** .5


def massHorizFrame(w):
    eccentricity = (1 - w ** 2 / length ** 2) ** .5
    circumference = 2 * length * \
        quad(lambda theta: (1 - eccentricity ** 2 *
                            np.sin(theta) ** 2) ** .5, 0, np.pi / 2)[0]
    return circumference * supportWidth * thickness * teakDensity


def massBottom(w):
    '''Add up the supports in the bottom to find the mass.'''
    midSkel = w * np.pi * ellipse(w, 0)
    midSkelMass = midSkel * supportWidth * thickness * teakDensity
    quarterSkel = np.pi * (ellipse(w, length / 2))
    quarterSkelMass = quarterSkel * supportWidth * thickness * teakDensity
    return quarterSkelMass * 2 + midSkelMass + massHorizFrame(w)

def circCg(w, r):
    cgOuter = (w - r / 3)
    mOuter = np.pi * 


def bottomCg(w):
    '''returns the center of mass for bottom of boat'''
    k
