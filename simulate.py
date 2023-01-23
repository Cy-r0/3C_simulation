#!/usr/bin/env python3
"""Chromosome Conformation Capture (3C) simulation"""


__author__ = 'Ciro Cursio'


import math

import numpy as np
import matplotlib.pyplot as plt


def generate_filament(angle_sigma=0.2):
    """Generate a sequence of 2D points with these constraints:
    1. Each pair of consecutive points is connected (i.e. there's a small
        constant distance between them);
    2. The direction of the vector between consecutive points is mostly random,
        but not completely: to encourage smooth curves, the direction between
        each pair is highly correlated with the direction between neighbouring
        pairs.
    Args:
        angle_sigma: variance of change in angle, expressed as fraction of 2pi
            (higher sigma = angle changes more sharply).
    """
    l = 1000
    distance = 1
    filament = np.zeros((l, 2))
    angle = 2 * math.pi * np.random.random()
    for i in range(1, l):
        v = [math.cos(angle) * distance, math.sin(angle) * distance]
        filament[i] = filament[i-1] + v
        angle += (np.random.random() - 0.5) * 2 * math.pi * angle_sigma
    return filament


def show_filament(filament):
    plt.scatter(filament[:, 0], filament[:, 1], marker='.')
    plt.show()


if __name__ == '__main__':
    f = generate_filament()
    show_filament(f)
