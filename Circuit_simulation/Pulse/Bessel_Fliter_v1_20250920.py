import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

import schemdraw
import schemdraw.elements as elm

with schemdraw.Drawing() as d:
    d.config(unit=2)
    d += (vin := elm.SourceV().label('Vin\n(50Ω)'))
    d += elm.Line().right()
    d.push()
    d += elm.Capacitor().down().label('C1=578pF')
    d += elm.Ground()
    d.pop()
    d += elm.Inductor().right().label('L2=1.69µH')
    d.push()
    d += elm.Capacitor().down().label('C3=719pF')
    d += elm.Ground()
    d.pop()
    d += elm.Inductor().right().label('L4=1.69µH')
    d.push()
    d += elm.Capacitor().down().label('C5=578pF')
    d += elm.Ground()
    d.pop()
    d += elm.Resistor().right().label('R_load\n(50Ω)')
    d += elm.SourceV().right().reverse().label('Vout').open()
    d.draw()