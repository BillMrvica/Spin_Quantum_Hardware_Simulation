# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 08:57:54 2023

@author: vjohn
"""

import numpy as np


class QuantumDot:
    # Constants
    MU_B = 9.274009994e-24  # Bohr magneton in J/T
    PLANCK = 6.62607015e-34  # Planck constant in J/Hz

    def __init__(self, g_factor, phi, theta):
        self.g_factor = g_factor
        self.phi = phi
        self.theta = theta
        self.shuttle_wait = 10
        self.shuttle_sweep = np.linspace(0, 500, 101)

    @property
    def axis(self):
        # Convert spherical coordinates to cartesian coordinates
        x = np.sin(self.theta) * np.cos(self.phi)
        y = np.sin(self.theta) * np.sin(self.phi)
        z = np.cos(self.theta)
        return x, y, z

    def lamor_freq(self, B):
        # Compute Larmor frequency using the provided formula
        return (self.g_factor * QuantumDot.MU_B * B) / QuantumDot.PLANCK

    def lamor_period(self, B):
        return 1/self.lamor_freq(B)
