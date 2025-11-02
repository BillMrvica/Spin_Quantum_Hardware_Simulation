# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 09:02:13 2023

@author: vjohn
"""

from array_343.figures.hopping_spins_paper.figureS24.QuantumDot import QuantumDot

#%%
class QuantumDotArray:
    def __init__(self, num_dots):
        self._num_dots = num_dots
        for i in range(num_dots):
            setattr(self, f"Q{i+1}", QuantumDot(None,
                                                None,
                                                None))

    def g_factors(self, dots=None):
        if dots is None:
            dots = range(1, self._num_dots+1)
        return [getattr(self, f"Q{i}").g_factor for i in dots]

    def phis(self, dots=None):
        if dots is None:
            dots = range(1, self._num_dots+1)
        return [getattr(self, f"Q{i}").phi for i in dots]

    def thetas(self, dots=None):
        if dots is None:
            dots = range(1, self._num_dots+1)
        return [getattr(self, f"Q{i}").theta for i in dots]

    def axes(self, dots=None):
        if dots is None:
            dots = range(1, self._num_dots+1)
        return [getattr(self, f"Q{i}").axis for i in dots]

    # Implementing the function to display all attributes in a compact and easy-to-read way
    def print_readable_snapshot(self, skip=['shuttle_sweep']):
        for i in range(self._num_dots):
            dot = getattr(self, f"Q{i+1}")
            print(f"--- Quantum Dot Q{i+1} ---")
            for key, value in vars(dot).items():
                if key not in skip:
                    print(f"{key}: {value}")
            print('_______________________')
