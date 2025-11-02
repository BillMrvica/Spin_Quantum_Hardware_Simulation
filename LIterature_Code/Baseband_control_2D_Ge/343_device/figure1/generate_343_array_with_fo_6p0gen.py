# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 15:47:21 2022

@author: valentinjohn
"""

# %% path management

import sys
import os

main_file = sys.modules['__main__'].__file__
script_dir = os.path.dirname(os.path.abspath(main_file))

if script_dir not in sys.path:
    sys.path.append(script_dir)

path = script_dir

# %% Imports
import gdstk
import numpy as np
from copy import copy
import math
from device_generator_343 import (gen_poly, gen_sensor, centroid,
                                  fanout_generator)

# %% generation and version
generation = 6
version = 0

# %% main array
file_name = f'343_array_{generation}p{version}gen.gds'

# %%% Input parameters

n_col = 3
n_row = 2

# Sizing all in um
spacing = 0.195 # dot-to-dot distance
plunger_width = 0.13
bar_width = 0.05 # barrier width
bar_width_sens = 0.06
bar_length = 0.074
bar_length_sens = 0.07
bar_offset = -0.002
chip_width = 4000 # full width of chip
margin = 300 # TODO
bondpad_spacing = 100 # distance between bondpads
bondpad_length = 400 # length of a bondpad
pad_int_fo_length = 0.1 # pad seperating small features below 0.1 um from larger fanout

min_feat_size = 0.03
gap_barrier = 0
gap_screening = 0.01
bend_radius = 0 #0.05

x_max = 0.95
y_max = 0.7
x_max2 = x_max + 1.3
y_max2 = y_max + 1.3
con_pl_width = 0.04
con_bar_width = 0.04
screen_width = con_bar_width + 0.01
sd_width = 0.05

# CAD layer numbers
pl_layer = 21
vbar_layer = 5
hbar_layer = 31


wbar_layer = 51
bbar_layer = 61
ohm_layer = 3

# conversion between small and large feature layer numbers
small_big_layer_dict = {3: 4, 5: 6, 21: 20, 31: 30, 51: 50, 61: 60}

# Calculate distances
cell_sp = n_col*spacing
use_width = chip_width - 2 * margin
vbar_length = n_row*spacing
hbar_length = n_col*spacing
spacing_hv = 2**0.5*spacing
feature_gap = (spacing - plunger_width - bar_width)/2
plunger_gap_long = spacing_hv - plunger_width
plunger_gap_short = spacing - plunger_width
con_bar_pos = (spacing_hv-plunger_width+(screen_width))/4

# Initialize gdstk library
lib = gdstk.Library()

# The design consists out of two sublattices, to acommodate the shifted plunger
# cells and alternating barrier order. A proper array has the form ABABA etc.
# (A and B the two sublattices)

# %%% Plungers
layer_b1 = hbar_layer
layer_b2 = vbar_layer
layer_b3 = wbar_layer
pl_layer = pl_layer

pl_cell = gdstk.Cell('plunger_v1') # define cell
plunger_points = gen_poly(8) # define polygon
plunger = gdstk.Polygon(plunger_points, layer = pl_layer) # add poly to cell
pl_cell.add(plunger)
plunger.scale(0.5/np.cos(np.pi/8)*plunger_width) # scale plunger to right size


overlap = (plunger_width + bar_width - spacing)/2

unit_cell = gdstk.Cell("unit_cell")
plunger_cell = gdstk.Cell("plunger_cell")

plunger_cell.add(gdstk.Reference(pl_cell, (-(n_col-1)/2*spacing_hv, -(n_row-1)/2*spacing_hv),
                              columns = n_col,
                              rows = n_row,
                              spacing=(spacing_hv, spacing_hv)
                              ),
              gdstk.Reference(pl_cell, (-(n_col)/2*spacing_hv, -(n_row-2)/2*spacing_hv),
                              columns = n_col+1,
                              rows = n_row-1,
                              spacing=(spacing_hv, spacing_hv)
                              )
              )

plunger_cell.flatten()

unit_cell.add(gdstk.Reference(plunger_cell))
unit_cell.flatten()

# %%%Inside barriers
bar_in_pos_cell = gdstk.Cell('bar_in_pos') # define cell
bar_in_neg_cell = gdstk.Cell('bar_in_neg') # define cell


for i in range(n_col):
    (x,y) = ((i-3/4)*spacing_hv+bar_offset, spacing_hv/4+bar_offset)

    bar_in_pos = gdstk.Polygon([(-bar_length/2, -bar_width/2),
                                (bar_length/2, -bar_width/2),
                                (bar_length/2+bar_width/4, -bar_width/8),
                                (bar_length/2+bar_width/4, bar_width/8),
                                (bar_length/2, bar_width/2),
                                (-bar_length/2, bar_width/2),
                                (-bar_length/2-bar_width/2, bar_width/8),
                                (-bar_length/2-bar_width/2, -bar_width/8)],
                               layer = layer_b3)
    bar_in_pos.set_property('name', 'bar_in')
    bar_in_pos.rotate(np.pi/4)
    bar_in_pos.translate(x, y)
    bar_in_pos_cell.add(bar_in_pos)

    bar_in_neg = gdstk.Polygon([(-bar_length/2, -bar_width/2),
                                (bar_length/2, -bar_width/2),
                                (bar_length/2+bar_width/2, -bar_width/8),
                                (bar_length/2+bar_width/2, bar_width/8),
                                (bar_length/2, bar_width/2),
                                (-bar_length/2, bar_width/2),
                                (-bar_length/2-bar_width/4, bar_width/8),
                                (-bar_length/2-bar_width/4, -bar_width/8),],
                               layer = layer_b2)
    bar_in_neg.set_property('name', 'bar_in')
    bar_in_neg.rotate(-np.pi/4)
    bar_in_neg.translate(-x, y)
    bar_in_neg_cell.add(bar_in_neg)

    bar_in_pos2 = bar_in_pos.copy()
    bar_in_pos2.mirror((1,0))
    bar_in_pos2.mirror((0,1))
    bar_in_pos_cell.add(bar_in_pos2)

    bar_in_neg2 = bar_in_neg.copy()
    bar_in_neg2.mirror((1,0))
    bar_in_neg2.mirror((0,1))
    bar_in_neg_cell.add(bar_in_neg2)


# %%% Outside barriers
bar_in_hor_cell = gdstk.Cell('bar_in_hor') # define cell
bar_in_ver_cell = gdstk.Cell('bar_in_ver') # define cell

bar_in_hor_top = gdstk.rectangle((-bar_length_sens/2, -bar_width_sens/2), (bar_length_sens/2, bar_width_sens/2), layer = layer_b1)
bar_in_hor_top.translate(0, spacing_hv/2+plunger_width/2+bar_width_sens/2+feature_gap)
bar_in_hor_cell.add(bar_in_hor_top)

bar_in_hor_bot = bar_in_hor_top.copy()
bar_in_hor_bot.mirror((1,0))
bar_in_hor_cell.add(bar_in_hor_bot)

bar_in_ver_left = gdstk.rectangle((bar_width_sens/2, bar_length_sens/2), (-bar_width_sens/2, -bar_length_sens/2), layer = layer_b1)
bar_in_ver_left.translate(spacing_hv*3/2+(plunger_width+bar_width_sens)/2+feature_gap, 0)
bar_in_ver_cell.add(bar_in_ver_left)

bar_in_ver_right = bar_in_ver_left.copy()
bar_in_ver_right.mirror((0,1))
bar_in_ver_cell.add(bar_in_ver_right)

unit_cell.add(gdstk.Reference(bar_in_pos_cell),
              gdstk.Reference(bar_in_neg_cell),
              gdstk.Reference(bar_in_hor_cell),
              gdstk.Reference(bar_in_ver_cell)
              )

# %%% Device consisting of unit cell

device = gdstk.Cell("device")
device.add(gdstk.cross((0,0), 3*plunger_width/4, 0.01),
           gdstk.Reference(unit_cell)
          )

device.flatten()

# %%% Sensor
sensor_width = 0.16
sensor_scale = (0.8,1.25)

sensor_position = ((0, spacing_hv/2+(sensor_scale[1]*sensor_width+plunger_width)/2+bar_width_sens+2*feature_gap),
                    (1.5*spacing_hv+(sensor_width+plunger_width)/2+bar_width_sens+2*feature_gap, 0),
                   (0, -(spacing_hv/2+(sensor_scale[1]*sensor_width+plunger_width)/2+bar_width_sens+2*feature_gap)),
                   (-(1.5*spacing_hv+(sensor_width+plunger_width)/2+bar_width_sens+2*feature_gap), 0))

sensor_width = 0.160
barrier_spacing = 0.008
ring_barrier_width = 0.03
bar_sens_width = 0.04
n_sharp = 0

sensor_north, pos_sd_north = gen_sensor(orientation_sensor='N',
                                        orientation_source='W',
                                        orientation_drain='E',
                                        sensor_width=sensor_width, bar_width=bar_sens_width,
                                        bar_length=bar_length_sens, sensor_sd_width = 0,
                                        layer_bar_sen=(layer_b3, layer_b2),
                                        bar_sou_position_offset=(-feature_gap,0.03),
                                        bar_dra_position_offset=(feature_gap,0.05),
                                        bar_sharp_source = bar_sens_width,
                                        bar_sharp_drain = bar_sens_width,
                                        # scale = sensor_scale
                                        )
sensor_east, pos_sd_east = gen_sensor(orientation_sensor='E',
                                      orientation_source='N',
                                      orientation_drain='S',
                                      sensor_width=sensor_width, bar_width=bar_sens_width,
                                      bar_length=bar_length_sens, sensor_sd_width = 0,
                                      layer_bar_sen=(layer_b2, layer_b3),
                                      source_position_offset=(0.01,0.),
                                      drain_position_offset=(0.04,0.0),
                                      bar_sou_position_offset=(-0.01, feature_gap),
                                      bar_dra_position_offset=(0.035, -feature_gap),
                                      bar_sharp_drain = 0.03
                                      )
sensor_south, pos_sd_south = gen_sensor(orientation_sensor='S',
                                        orientation_source='E',
                                        orientation_drain='W',
                                        sensor_width=sensor_width, bar_width=bar_sens_width,
                                        bar_length=bar_length_sens, sensor_sd_width = 0,
                                        layer_bar_sen=(layer_b3, layer_b2),
                                        bar_dra_position_offset=(-feature_gap,-0.05),
                                        bar_sou_position_offset=(feature_gap,-0.03),
                                        bar_sharp_source = bar_sens_width,
                                        bar_sharp_drain = bar_sens_width,
                                        )
sensor_west, pos_sd_west = gen_sensor(orientation_sensor='W',
                                      orientation_source='S',
                                      orientation_drain='N',
                                      sensor_width=sensor_width, bar_width=bar_sens_width,
                                      bar_length=bar_length_sens, sensor_sd_width = 0,
                                      layer_bar_sen=(layer_b2, layer_b3),
                                      drain_position_offset=(-0.04,0.0),
                                      source_position_offset=(-0.01,0.0),
                                      bar_dra_position_offset=(-0.035,+feature_gap),
                                      bar_sou_position_offset=(0, -feature_gap),
                                      bar_sharp_drain = 0.03
                                      )

sensor_cell = gdstk.Cell('sensor')
sensor_cell.add(gdstk.Reference(sensor_north, sensor_position[0]))
sensor_cell.add(gdstk.Reference(sensor_east, sensor_position[1]))
sensor_cell.add(gdstk.Reference(sensor_south, sensor_position[2]))
sensor_cell.add(gdstk.Reference(sensor_west, sensor_position[3]))

sensor_cell.flatten()

# %%% Connection lines positions
n_con_dict = {'north':14,
              'east':13,
              'south':14,
              'west':13}
con_pos_dict = {'north':[],
              'east':[],
              'south':[],
              'west':[]}

edge_space = 0.06
for x in np.linspace(-(x_max-edge_space),(x_max-edge_space),n_con_dict['north']):
    pos = (x, y_max)
    con_pos_dict['north'].append(pos)

for x in np.linspace(-(x_max-edge_space),(x_max-edge_space),n_con_dict['south']):
    pos = (x, -y_max)
    con_pos_dict['south'].append(pos)

for y in np.linspace(-(y_max-edge_space),(y_max-edge_space),n_con_dict['east']):
    pos = (x_max, y)
    con_pos_dict['east'].append(pos)

for y in np.linspace(-(y_max-edge_space),(y_max-edge_space),n_con_dict['west']):
    pos = (-x_max, y)
    con_pos_dict['west'].append(pos)

con_pos_dict_2 = {'north':[],
              'east':[],
              'south':[],
              'west':[]}

edge_space2 = 0.3
for x in np.linspace(-(x_max2-edge_space2),(x_max2-edge_space2),n_con_dict['north']):
    pos = (x, y_max2)
    con_pos_dict_2['north'].append(pos)

for x in np.linspace(-(x_max2-edge_space2),(x_max2-edge_space2),n_con_dict['south']):
    pos = (x, -y_max2)
    con_pos_dict_2['south'].append(pos)

for y in np.linspace(-(y_max2-edge_space2),(y_max2-edge_space2),n_con_dict['east']):
    pos = (x_max2, y)
    con_pos_dict_2['east'].append(pos)

for y in np.linspace(-(y_max2-edge_space2),(y_max2-edge_space2),n_con_dict['west']):
    pos = (-x_max2, y)
    con_pos_dict_2['west'].append(pos)

# %%% Connecting sensor
space_bar_sd = feature_gap


bar_con_cell = gdstk.Cell('bar_con')

# North sensor
con_pl_sens_north = gdstk.FlexPath([sensor_position[0],
                                    (sensor_position[0][0], sensor_position[0][1]+0.13),
                                    (sensor_position[0][0], sensor_position[0][1]+0.325),
                                    ],
                                    con_bar_width, simple_path=False, layer=pl_layer)
con_pl_sens_north.segment((con_pos_dict_2['north'][6][0], con_pos_dict_2['north'][6][1]-pad_int_fo_length), 0.1)
con_pl_sens_north.segment(con_pos_dict_2['north'][6], 0.1)
bar_con_cell.add(con_pl_sens_north)

source_north = gdstk.Polygon([(sensor_cell.polygons[3].points[0][0]-space_bar_sd, sensor_cell.polygons[3].points[0][1]),
                                         (sensor_cell.polygons[3].points[1][0]-space_bar_sd, sensor_cell.polygons[3].points[1][1]),
                                         (sensor_cell.polygons[3].points[1][0]-space_bar_sd-0.01, sensor_cell.polygons[3].points[1][1]),
                                         (sensor_cell.polygons[3].points[1][0]-space_bar_sd-0.01-2**0.5/2*sd_width, sensor_cell.polygons[3].points[1][1]-2**0.5/2*sd_width)
                                         ], layer = ohm_layer)
bar_con_cell.add(source_north)

con_sd_sens_north_sou = gdstk.FlexPath([centroid(source_north.points[2:4]),
                                        (centroid(source_north.points[2:4])[0]-0.04,centroid(source_north.points[2:4])[1]+0.04),
                                        con_pos_dict['north'][4]],
                                       sd_width, simple_path=False, layer=ohm_layer)
con_sd_sens_north_sou.segment((con_pos_dict_2['north'][4][0], con_pos_dict_2['north'][4][1]-pad_int_fo_length), 0.1)
con_sd_sens_north_sou.segment(con_pos_dict_2['north'][4], 0.1)
bar_con_cell.add(con_sd_sens_north_sou)

drain_north = gdstk.Polygon([(sensor_cell.polygons[4].points[3][0]+space_bar_sd, sensor_cell.polygons[4].points[3][1]+0.005),
                                         (sensor_cell.polygons[4].points[2][0]+space_bar_sd, sensor_cell.polygons[4].points[2][1]+0.03),
                                         (sensor_cell.polygons[4].points[2][0]+space_bar_sd+0.01, sensor_cell.polygons[4].points[2][1]+0.03),
                                         (sensor_cell.polygons[4].points[2][0]+space_bar_sd+0.01+2**0.5/2*sd_width, sensor_cell.polygons[4].points[2][1]-2**0.5/2*sd_width+0.03)
                                         ], layer = ohm_layer)
bar_con_cell.add(drain_north)

con_sd_sens_north_dra = gdstk.FlexPath([centroid(drain_north.points[2:4]),
                                        (centroid(drain_north.points[2:4])[0]+0.02,centroid(drain_north.points[2:4])[1]+0.02),
                                        (con_pos_dict['north'][-5][0]-0.08, con_pos_dict['north'][-5][1])],
                                       sd_width, simple_path=False, layer=ohm_layer)
con_sd_sens_north_dra.segment((con_pos_dict_2['north'][-5][0], con_pos_dict_2['north'][-5][1]-pad_int_fo_length), 0.1)
con_sd_sens_north_dra.segment(con_pos_dict_2['north'][-5], 0.1)
bar_con_cell.add(con_sd_sens_north_dra)


con_bar_sens_north_sou = gdstk.FlexPath([centroid(sensor_cell.polygons[3].points),
                                         (centroid(sensor_cell.polygons[3].points)[0],
                                          centroid(sensor_cell.polygons[3].points)[1]+0.1),
                                        ],
                                       con_bar_width, simple_path=False, layer=layer_b3)
con_bar_sens_north_sou.segment((con_pos_dict_2['north'][5][0], con_pos_dict_2['north'][5][1]-pad_int_fo_length), 0.1)
con_bar_sens_north_sou.segment(con_pos_dict_2['north'][5], 0.1)
con_bar_sens_north_dra = gdstk.FlexPath([centroid(sensor_cell.polygons[4].points),
                                         (centroid(sensor_cell.polygons[4].points)[0],
                                          centroid(sensor_cell.polygons[4].points)[1]+4*con_bar_width)
                                         ],
                                       con_bar_width, simple_path=False, layer=layer_b2)
con_bar_sens_north_dra.segment((con_pos_dict_2['north'][8][0], con_pos_dict_2['north'][8][1]-pad_int_fo_length), 0.1)
con_bar_sens_north_dra.segment(con_pos_dict_2['north'][8], 0.1)
bar_con_cell.add(con_bar_sens_north_sou,
                 con_bar_sens_north_dra)

# East sensor
con_pl_sens_east = gdstk.FlexPath([sensor_position[1],
                                   (sensor_position[1][0]+0.15, sensor_position[1][1]),
                                   (sensor_position[1][0]+0.3, sensor_position[1][1]),
                                   (sensor_position[1][0]+0.33, sensor_position[1][1])
                                    ],
                                    con_bar_width, simple_path=False, layer=pl_layer)
con_pl_sens_east.segment((con_pos_dict_2['east'][6][0]-pad_int_fo_length, con_pos_dict_2['east'][6][1]), 0.1)
con_pl_sens_east.segment(con_pos_dict_2['east'][6], 0.1)
bar_con_cell.add(con_pl_sens_east)

drain_east = gdstk.Polygon([(sensor_cell.polygons[18+n_sharp].points[0][0]+0.01, sensor_cell.polygons[18+n_sharp].points[0][1]+space_bar_sd),
                                         (sensor_cell.polygons[18+n_sharp].points[1][0]+0.01, sensor_cell.polygons[18+n_sharp].points[1][1]+space_bar_sd),
                                         (sensor_cell.polygons[18+n_sharp].points[1][0]+0.01, sensor_cell.polygons[18+n_sharp].points[1][1]+space_bar_sd+0.01),
                                         (sensor_cell.polygons[18+n_sharp].points[1][0]+0.01-2**0.5/2*sd_width, sensor_cell.polygons[18+n_sharp].points[1][1]+space_bar_sd+0.01+2**0.5/2*sd_width)
                                         ], layer = ohm_layer)
bar_con_cell.add(drain_east)

source_east = gdstk.Polygon([(sensor_cell.polygons[19+n_sharp].points[3][0], sensor_cell.polygons[19+n_sharp].points[3][1]-space_bar_sd),
                                         (sensor_cell.polygons[19+n_sharp].points[2][0], sensor_cell.polygons[19+n_sharp].points[2][1]-space_bar_sd),
                                         (sensor_cell.polygons[19+n_sharp].points[2][0], sensor_cell.polygons[19+n_sharp].points[2][1]-space_bar_sd-0.01),
                                         (sensor_cell.polygons[19+n_sharp].points[2][0]-2**0.5/2*sd_width, sensor_cell.polygons[19+n_sharp].points[2][1]-space_bar_sd-0.01-2**0.5/2*sd_width)
                                         ], layer = ohm_layer)
bar_con_cell.add(source_east)

con_sd_sens_east_sou = gdstk.FlexPath([centroid(source_east.points[2:4]),
                                       (centroid(source_east.points[2:4])[0]+0.01,centroid(source_east.points[2:4])[1]-0.01),
                                       con_pos_dict['east'][3]],
                                       sd_width, simple_path=False, layer=ohm_layer)
con_sd_sens_east_sou.segment((con_pos_dict_2['east'][3][0]-pad_int_fo_length, con_pos_dict_2['east'][3][1]), 0.1)
con_sd_sens_east_sou.segment(con_pos_dict_2['east'][3], 0.1)
bar_con_cell.add(con_sd_sens_east_sou)

con_sd_sens_east_dra = gdstk.FlexPath([centroid(drain_east.points[2:4]),
                                       (centroid(drain_east.points[2:4])[0]+0.02,centroid(drain_east.points[2:4])[1]+0.02),
                                       (con_pos_dict['east'][-5][0], con_pos_dict['east'][-5][1])],
                                       sd_width, simple_path=False, layer=ohm_layer)
con_sd_sens_east_dra.segment((con_pos_dict_2['east'][-5][0]-pad_int_fo_length, con_pos_dict_2['east'][-5][1]), 0.1)
con_sd_sens_east_dra.segment(con_pos_dict_2['east'][-5], 0.1)
bar_con_cell.add(con_sd_sens_east_dra)


con_bar_sens_east_sou = gdstk.FlexPath([centroid(sensor_cell.polygons[19+n_sharp].points),
                                        (centroid(sensor_cell.polygons[19+n_sharp].points)[0]+2*bar_length_sens/2,
                                         centroid(sensor_cell.polygons[19+n_sharp].points)[1]),
                                        con_pos_dict['east'][4]],
                                       con_bar_width, simple_path=False, layer=layer_b3)
con_bar_sens_east_sou.segment((con_pos_dict_2['east'][4][0]-pad_int_fo_length, con_pos_dict_2['east'][4][1]), 0.1)
con_bar_sens_east_sou.segment(con_pos_dict_2['east'][4], 0.1)

con_bar_sens_east_dra = gdstk.FlexPath([centroid(sensor_cell.polygons[18+n_sharp].points),
                                        (centroid(sensor_cell.polygons[18+n_sharp].points)[0]+0.3,
                                         centroid(sensor_cell.polygons[18+n_sharp].points)[1]),
                                        ],
                                       con_bar_width, simple_path=False, layer=layer_b2)
con_bar_sens_east_dra.segment((con_pos_dict_2['east'][7][0]-pad_int_fo_length, con_pos_dict_2['east'][7][1]), 0.1)
con_bar_sens_east_dra.segment(con_pos_dict_2['east'][7], 0.1)

bar_con_cell.add(con_bar_sens_east_sou,
                 con_bar_sens_east_dra)

# South sensor
drain_south = drain_north.copy()
drain_south.mirror((1,0))
drain_south.mirror((0,1))
source_south = source_north.copy()
source_south.mirror((1,0))
source_south.mirror((0,1))

con_pl_sens_south = con_pl_sens_north.copy()
con_pl_sens_south.mirror((1,0))
con_pl_sens_south.mirror((0,1))
con_sd_sens_south_sou = con_sd_sens_north_sou.copy()
con_sd_sens_south_sou.mirror((1,0))
con_sd_sens_south_sou.mirror((0,1))
con_sd_sens_south_dra = con_sd_sens_north_dra.copy()
con_sd_sens_south_dra.mirror((1,0))
con_sd_sens_south_dra.mirror((0,1))
con_bar_sens_south_sou = con_bar_sens_north_sou.copy()
con_bar_sens_south_sou.mirror((1,0))
con_bar_sens_south_sou.mirror((0,1))
con_bar_sens_south_dra = con_bar_sens_north_dra.copy()
con_bar_sens_south_dra.mirror((1,0))
con_bar_sens_south_dra.mirror((0,1))

bar_con_cell.add(con_pl_sens_south,
                 con_sd_sens_south_sou,
                 con_sd_sens_south_dra,
                 con_bar_sens_south_sou,
                 con_bar_sens_south_dra,
                 drain_south,
                 source_south)

# West sensor
drain_west = drain_east.copy()
drain_west.mirror((1,0))
drain_west.mirror((0,1))
source_west = source_east.copy()
source_west.mirror((1,0))
source_west.mirror((0,1))

con_pl_sens_west = con_pl_sens_east.copy()
con_pl_sens_west.mirror((0,1))
con_pl_sens_west.mirror((1,0))
con_sd_sens_west_sou = con_sd_sens_east_sou.copy()
con_sd_sens_west_sou.mirror((0,1))
con_sd_sens_west_sou.mirror((1,0))
con_sd_sens_west_dra = con_sd_sens_east_dra.copy()
con_sd_sens_west_dra.mirror((0,1))
con_sd_sens_west_dra.mirror((1,0))
con_bar_sens_west_sou = con_bar_sens_east_sou.copy()
con_bar_sens_west_sou.mirror((0,1))
con_bar_sens_west_sou.mirror((1,0))
con_bar_sens_west_dra = con_bar_sens_east_dra.copy()
con_bar_sens_west_dra.mirror((0,1))
con_bar_sens_west_dra.mirror((1,0))


bar_con_cell.add(con_pl_sens_west,
                 con_sd_sens_west_sou,
                 con_sd_sens_west_dra,
                 con_bar_sens_west_sou,
                 con_bar_sens_west_dra,
                 drain_west,
                 source_west)

# %%% Connecting plungers
# Plunger 11
con_pl_11 = gdstk.FlexPath([(-spacing_hv,spacing_hv/2),
                            (-spacing_hv-0.1,spacing_hv/2+0.07),
                            (-spacing_hv-0.15,spacing_hv/2+0.10),
                            (-7/4*spacing_hv-0.01, 4/4*spacing_hv+0.017),
                            con_pos_dict['west'][-1]
                            ],
                           con_bar_width, simple_path=False, layer=pl_layer, bend_radius=bend_radius)
con_pl_11.segment((con_pos_dict_2['west'][-1][0]+pad_int_fo_length, con_pos_dict_2['west'][-1][1]), 0.1)
con_pl_11.segment(con_pos_dict_2['west'][-1], 0.1)
bar_con_cell.add(con_pl_11)
# Plunger 12
con_pl_12 = gdstk.FlexPath([(0,spacing_hv/2),
                            (-spacing_hv/2,spacing_hv/2),
                            (-spacing_hv/2,spacing_hv/2+plunger_width/2-0.006),
                            (-spacing_hv*3/4-0.02,spacing_hv*3/4+plunger_width/2-0.006+0.02),
                            (-spacing_hv*3/4-0.06,spacing_hv*3/4+plunger_width/2+0.035),
                            (-spacing_hv*5/4-0.04,spacing_hv*5/4+plunger_width/2+0.0),
                            con_pos_dict['north'][1]],
                           con_bar_width, simple_path=False, layer=pl_layer, bend_radius=0)
con_pl_12.segment((con_pos_dict_2['north'][1][0], con_pos_dict_2['north'][1][1]-pad_int_fo_length), 0.1)
con_pl_12.segment((con_pos_dict_2['north'][1][0], con_pos_dict_2['north'][1][1]), 0.1)
bar_con_cell.add(con_pl_12)
# Plunger 13
con_pl_13 = gdstk.FlexPath([(spacing_hv,spacing_hv/2),
                            (5.5/4*spacing_hv, 3.5/4*spacing_hv),
                            (6/4*spacing_hv+0.07, 5/4*spacing_hv),
                            (8/4*spacing_hv+0.05, 7/4*spacing_hv-0.03),
                            (con_pos_dict['east'][-1][0], con_pos_dict['east'][-1][1]+0.04)],
                           con_bar_width, simple_path=False, layer=pl_layer, bend_radius=bend_radius)
con_pl_13.segment((con_pos_dict_2['east'][-1][0]-pad_int_fo_length, con_pos_dict_2['east'][-1][1]), 0.1)
con_pl_13.segment(con_pos_dict_2['east'][-1], 0.1)
bar_con_cell.add(con_pl_13)
# Plunger 20
# Plunger 11
con_pl_20 = gdstk.FlexPath([(-1.5*spacing_hv-0.025, 0),
                            (-1.5*spacing_hv-0.025, -spacing_hv/2+0.04),
                            (-8/4*spacing_hv-0.00, -3/4*spacing_hv-con_bar_width*3/4+0.01),
                            (-8/4*spacing_hv-0.05, -3/4*spacing_hv-con_bar_width*3/4-0.03),
                            (-9/4*spacing_hv-0.11, -4/4*spacing_hv-con_bar_width*3/4-0.05),
                            con_pos_dict['west'][2]],
                           con_bar_width, simple_path=False, layer=pl_layer, bend_radius=bend_radius)
con_pl_20.segment((con_pos_dict_2['west'][2][0]+pad_int_fo_length, con_pos_dict_2['west'][2][1]), 0.1)
con_pl_20.segment(con_pos_dict_2['west'][2], 0.1)
bar_con_cell.add(con_pl_20)
# Plunger 23
con_pl_23 = con_pl_20.copy()
con_pl_23.mirror((0,1))
con_pl_23.mirror((1,0))
bar_con_cell.add(con_pl_23)
# Plunger 22
con_pl_22 = gdstk.FlexPath([(spacing_hv/2,0),
                            (spacing_hv/2,spacing_hv/2+plunger_width/2-0.006),
                            (spacing_hv*3/4+0.06,spacing_hv*3/4+plunger_width/2-0.006+0.06),
                            (spacing_hv*3/4+0.15,spacing_hv*3/4+plunger_width/2+0.135),
                            (spacing_hv*5/4+0.15,spacing_hv*5/4+plunger_width/2+0.135),
                            con_pos_dict['north'][-3]],
                           con_bar_width, simple_path=False, layer=pl_layer, bend_radius=bend_radius)
con_pl_22.segment((con_pos_dict_2['north'][-2][0], con_pos_dict_2['north'][-2][1]-pad_int_fo_length), 0.1)
con_pl_22.segment(con_pos_dict_2['north'][-2], 0.1)
bar_con_cell.add(con_pl_22)

# Plunger 21
con_pl_21 = con_pl_22.copy()
con_pl_21.mirror((0,1))
con_pl_21.mirror((1,0))
bar_con_cell.add(con_pl_21)

# Plunger 31
con_pl_31 = con_pl_13.copy()
con_pl_31.mirror((1,0))
con_pl_31.mirror((0,1))
bar_con_cell.add(con_pl_31)
bar_con_cell.add(con_pl_31)
# Plunger 32
con_pl_32 = con_pl_12.copy()
con_pl_32.mirror((1,0))
con_pl_32.mirror((0,1))
bar_con_cell.add(con_pl_32)
# Plunger 33
con_pl_33 = con_pl_11.copy()
con_pl_33.mirror((1,0))
con_pl_33.mirror((0,1))
bar_con_cell.add(con_pl_33)

# %%% Connecting barriers
# sensor barrier north

con_bar_in_hor_north = gdstk.FlexPath([(0, spacing_hv-plunger_gap_long/2+(bar_width_sens+bar_width)/2-con_bar_width/2),
                                       (0.03, spacing_hv-plunger_gap_long/2+(bar_width_sens+bar_width)/2-con_bar_width/2),
                                       (0.035+1/4*spacing_hv, 6/4*spacing_hv-plunger_gap_long/2),
                                       (0.11+2/4*spacing_hv, 6/4*spacing_hv-plunger_gap_long/2+bar_width_sens+0.08),
                                       (con_pos_dict['north'][-4][0]-0.1, con_pos_dict['north'][-4][1])],
                                      con_bar_width, simple_path=False, layer=layer_b1, bend_radius=bend_radius)
con_bar_in_hor_north.segment((con_pos_dict_2['north'][-4][0], con_pos_dict_2['north'][-4][1]-pad_int_fo_length), 0.1)
con_bar_in_hor_north.segment(con_pos_dict_2['north'][-4], 0.1)
bar_con_cell.add(con_bar_in_hor_north)

# sensor barrier east
x1 = bar_in_ver_cell.polygons[0].points[1][0] + (bar_width-con_bar_width)/2
x2 = bar_in_ver_cell.polygons[0].points[2][0] - (bar_width-con_bar_width)/2
y1 = bar_in_ver_cell.polygons[0].points[1][1]
y2 = 0.05
con_bar_in_ver_east = gdstk.rectangle((x1,y1), (x2,y2), layer = bar_in_hor_cell.polygons[0].layer)

con_bar_in_hor_east = gdstk.FlexPath([(2*spacing_hv-plunger_gap_long/2+bar_width-con_bar_width/2, 0),
                                       (2*spacing_hv-plunger_gap_long/2+bar_width-con_bar_width/2, -0.07),
                                       (10/4*spacing_hv-plunger_gap_long/2+bar_width-con_bar_width/2, -(0.09+2/4*spacing_hv)),
                                       con_pos_dict['east'][2]],
                                      con_bar_width, simple_path=False, layer=layer_b1, bend_radius=bend_radius)
con_bar_in_hor_east.segment((con_pos_dict_2['east'][2][0]-pad_int_fo_length, con_pos_dict_2['east'][2][1]), 0.1)
con_bar_in_hor_east.segment(con_pos_dict_2['east'][2], 0.1)
bar_con_cell.add(con_bar_in_hor_east)

# sensor barrier south
con_bar_in_hor_south = con_bar_in_hor_north.copy()
con_bar_in_hor_south.mirror((0,1))
con_bar_in_hor_south.mirror((1,0))
bar_con_cell.add(con_bar_in_hor_south)

# sensor barrier west
con_bar_in_hor_west = con_bar_in_hor_east.copy()
con_bar_in_hor_west.mirror((0,1))
con_bar_in_hor_west.mirror((1,0))
bar_con_cell.add(con_bar_in_hor_west)

con_bar_11 = gdstk.FlexPath([(-5/4*spacing_hv, spacing_hv/4),
                             (-5/4*spacing_hv-1.7*con_bar_width/2, spacing_hv/4+1.7*con_bar_width/2),
                             (-6/4*spacing_hv-1.7*con_bar_width/2, 3/8*spacing_hv+1.7*con_bar_width/2),
                             (-7/4*spacing_hv-1.7*con_bar_width/2-0.03, 4/4*spacing_hv-0.03),
                             (-7/4*spacing_hv-1.7*con_bar_width/2-0.08, 4/4*spacing_hv+0.02),
                              con_pos_dict['west'][-2]],
                     con_bar_width, simple_path=False, layer=layer_b2, bend_radius=bend_radius)
con_bar_11.segment((con_pos_dict_2['west'][-2][0]+pad_int_fo_length, con_pos_dict_2['west'][-2][1]), 0.1)
con_bar_11.segment(con_pos_dict_2['west'][-2], 0.1)
bar_con_cell.add(con_bar_11)

con_bar_12 = gdstk.FlexPath([(-3/4*spacing_hv, spacing_hv/4),
                             (-1/2*spacing_hv-con_bar_pos, spacing_hv/4+con_bar_width/2+feature_gap/2),
                             (-1/2*spacing_hv-con_bar_pos, spacing_hv*3/4-con_bar_width/2-0.015+feature_gap/2),
                             (-1/2*spacing_hv-con_bar_pos-0.06, spacing_hv*3/4-con_bar_width/2-0.01+feature_gap/2+0.055),
                             (-1/2*spacing_hv-con_bar_pos-0.13, spacing_hv*3/4-con_bar_width/2-0.006+0.08),
                             con_pos_dict['north'][0]],
                             con_bar_width, simple_path=False, layer=layer_b3, bend_radius=0)
con_bar_12.segment((con_pos_dict_2['north'][0][0], con_pos_dict_2['north'][0][1]-pad_int_fo_length), 0.1)
con_bar_12.segment(con_pos_dict_2['north'][0], 0.1)
bar_con_cell.add(con_bar_12)

con_bar_13 = gdstk.FlexPath([(-1/4*spacing_hv, spacing_hv/4),
                             (-1/2*spacing_hv+con_bar_pos, spacing_hv/4+con_bar_width/2+feature_gap/2),
                             (-1/2*spacing_hv+con_bar_pos, 3/4*spacing_hv+0.013+feature_gap/2),
                             (-1/2*spacing_hv+con_bar_pos-0.147, 3/4*spacing_hv+0.013+feature_gap/2+0.147),
                             (-1/2*spacing_hv+con_bar_pos-0.2, 3/4*spacing_hv+0.27),
                             con_pos_dict['north'][3]],
                     con_bar_width, simple_path=False, layer=layer_b2, bend_radius=0)
con_bar_13.segment((con_pos_dict_2['north'][3][0], con_pos_dict_2['north'][3][1]-pad_int_fo_length), 0.1)
con_bar_13.segment(con_pos_dict_2['north'][3], 0.1)
bar_con_cell.add(con_bar_13)

con_bar_14 = gdstk.FlexPath([(1/4*spacing_hv, spacing_hv/4),
                             (1/2*spacing_hv-con_bar_pos, spacing_hv/4+con_bar_width/2+feature_gap/2),
                             (1/2*spacing_hv-con_bar_pos, 3/4*spacing_hv+feature_gap/2+0.022),
                             (3/4*spacing_hv+con_bar_width/2+feature_gap/2, 5/4*spacing_hv+feature_gap/2+0.022),
                             (5/4*spacing_hv+con_bar_width/2+feature_gap/2, 7/4*spacing_hv+feature_gap/2+0.022),
                             con_pos_dict['north'][-4]],
                     con_bar_width, simple_path=False, layer=layer_b3, bend_radius=bend_radius)
con_bar_14.segment((con_pos_dict_2['north'][-3][0], con_pos_dict_2['north'][-3][1]-pad_int_fo_length), 0.1)
con_bar_14.segment(con_pos_dict_2['north'][-3], 0.1)
bar_con_cell.add(con_bar_14)

con_bar_15 = gdstk.FlexPath([(3/4*spacing_hv, spacing_hv/4),
                             (1/2*spacing_hv+con_bar_pos, spacing_hv/4+con_bar_width/2+feature_gap/2),
                             (1/2*spacing_hv+con_bar_pos, spacing_hv*3/4-con_bar_width/2+feature_gap/2-0.015),
                             (1/2*spacing_hv+con_bar_pos+0.045, spacing_hv*3/4-con_bar_width/2+feature_gap/2-0.015+0.045),
                             (3/4*spacing_hv-con_bar_width/2-feature_gap/2+0.11, spacing_hv*3/4-con_bar_width/2+feature_gap/2-0.006+0.07),
                             (5/4*spacing_hv-con_bar_width/2-feature_gap/2+0.1, spacing_hv*5/4-con_bar_width/2+feature_gap/2-0.006+0.07),
                             (con_pos_dict['north'][-1][0]-0.1, con_pos_dict['north'][-1][1])],
                             con_bar_width, simple_path=False, layer=layer_b2, bend_radius=bend_radius)
con_bar_15.segment((con_pos_dict_2['north'][-1][0], con_pos_dict_2['north'][-1][1]-pad_int_fo_length), 0.1)
con_bar_15.segment(con_pos_dict_2['north'][-1], 0.1)
bar_con_cell.add(con_bar_15)

con_bar_16 = gdstk.FlexPath([(5/4*spacing_hv, spacing_hv/4),
                             (5/4*spacing_hv+1.7*con_bar_width/2, spacing_hv/4+1.7*con_bar_width/2),
                             (5/4*spacing_hv+1.7*con_bar_width/2+0.02, 2/4*spacing_hv+0.03),
                             (5/4*spacing_hv+1.7*con_bar_width/2+0.05, 2/4*spacing_hv+0.06),
                             (7/4*spacing_hv+1.7*con_bar_width/2+0.07, 4/4*spacing_hv+0.06),
                             (9/4*spacing_hv+1.7*con_bar_width/2+0.07, 5/4*spacing_hv+0.08),
                              con_pos_dict['east'][-2]],
                     con_bar_width, simple_path=False, layer=layer_b3, bend_radius=bend_radius)
con_bar_16.segment((con_pos_dict_2['east'][-2][0]-pad_int_fo_length, con_pos_dict_2['west'][-2][1]), 0.1)
con_bar_16.segment(con_pos_dict_2['east'][-2], 0.1)
bar_con_cell.add(con_bar_16)

con_bar_21 = con_bar_16.copy()
con_bar_21.mirror((1,0))
con_bar_21.mirror((0,1))
bar_con_cell.add(con_bar_21)

con_bar_22 = con_bar_15.copy()
con_bar_22 = con_bar_22.mirror((1,0))
con_bar_22 = con_bar_22.mirror((0,1))
bar_con_cell.add(con_bar_22)
bar_con_cell.add(con_bar_22)

con_bar_23 = con_bar_14.copy()
con_bar_23 = con_bar_23.mirror((1,0))
con_bar_23 = con_bar_23.mirror((0,1))
bar_con_cell.add(con_bar_23)

con_bar_24 = con_bar_13.copy()
con_bar_24 = con_bar_24.mirror((1,0))
con_bar_24 = con_bar_24.mirror((0,1))
bar_con_cell.add(con_bar_24)

con_bar_25 = con_bar_12.copy()
con_bar_25 = con_bar_25.mirror((1,0))
con_bar_25 = con_bar_25.mirror((0,1))
bar_con_cell.add(con_bar_25)

con_bar_26 = con_bar_11.copy()
con_bar_26.mirror((1,0))
con_bar_26.mirror((0,1))
bar_con_cell.add(con_bar_26)

# %%% Screening gates

screening_gate_cell = gdstk.Cell("screening_gate")
screening_gate_sensor_cell = gdstk.Cell("screening_gate_sensors")


con_pl_list = [con_pl_sens_north, con_pl_sens_east,
               con_pl_sens_south, con_pl_sens_west,
               con_pl_11, con_pl_12, con_pl_13,
               con_pl_20, con_pl_21, con_pl_22, con_pl_23,
               con_pl_31, con_pl_32, con_pl_33]

con_pl_screen_start = [1,1,
                       1,1,
                       1,1,1,
                       1,0,0,1,
                       1,1,1
                       ]

con_pl_screen_end = [3,4,
                     3,4,
                     4,6,3,
                     4,3,3,4,
                     3,6,4
                       ]


centre_screen_finger = {}
end_screen_finger = {}
screening_gate_dict = {}
n = 0
for con_pl in con_pl_list:
    if con_pl == con_pl_12 or con_pl == con_pl_32:
        ends = (0.03, 0.03)
    else: ends = (-0.02, 0.03)
    screening_gate_N_rect = gdstk.FlexPath(con_pl.path_spines()[0][con_pl_screen_start[n]:con_pl_screen_end[n]],
                                           screen_width, simple_path=False, layer=layer_b1, ends=ends)
    end_screen_finger[con_pl] = screening_gate_N_rect.path_spines()[0][-1]
    screening_gate_dict[con_pl] = screening_gate_N_rect
    n = n + 1

#  north-west screening gates
screening_pl_12 = gdstk.FlexPath(screening_gate_dict[con_pl_12].path_spines()[0]
                                   , screen_width, simple_path=False, layer=layer_b1, ends=(0.03,0))
screening_pl_12.segment((end_screen_finger[con_pl_12][0]+0.02,
                          end_screen_finger[con_pl_12][1]+0.05), con_bar_width)
screening_pl_12.segment(con_pos_dict['north'][2])

screening_pl_11 = gdstk.FlexPath(screening_gate_dict[con_pl_11].path_spines()[0][:-1]
                                   , screen_width+0.005, simple_path=False, layer=layer_b1, ends=(0,0))
screening_pl_11.segment(screening_gate_dict[con_pl_11].path_spines()[0][-1], screen_width+0.015)
screening_pl_11.segment(end_screen_finger[con_pl_12], screen_width++0.02)
screening_pl_11.segment((end_screen_finger[con_pl_12][0]+0.02,
                          end_screen_finger[con_pl_12][1]+0.05), con_bar_width)
screening_pl_11.segment(con_pos_dict['north'][2])
screening_pl_11.segment((con_pos_dict_2['north'][2][0], con_pos_dict_2['north'][2][1]-pad_int_fo_length), 0.1)
screening_pl_11.segment(con_pos_dict_2['north'][2], 0.1)

screening_pl_sens_north = gdstk.FlexPath(screening_gate_dict[con_pl_sens_north].path_spines()[0]
                                      , 0.06, simple_path=False, layer=layer_b1, ends=(-0.02,0))

screening_pl_sens_north.segment(end_screen_finger[con_pl_sens_north], screen_width)
screening_pl_sens_north.segment((end_screen_finger[con_pl_sens_north][0]+0.05, end_screen_finger[con_pl_sens_north][1]+0.02), con_bar_width+0.0)
screening_pl_sens_north.segment((con_pos_dict_2['north'][7][0], con_pos_dict_2['north'][7][1]-pad_int_fo_length), 0.1)
screening_pl_sens_north.segment(con_pos_dict_2['north'][7], 0.1)


#  north-east screening gates
screening_pl_22 = gdstk.FlexPath(screening_gate_dict[con_pl_22].path_spines()[0][:-1]
                                  , screen_width, simple_path=False, layer=layer_b1, ends=(-0.108,0))
screening_pl_22.segment(screening_gate_dict[con_pl_22].path_spines()[0][-1], screen_width)
screening_pl_22.segment(screening_gate_dict[con_pl_13].path_spines()[0], screen_width+0.015)
screening_pl_22.segment(end_screen_finger[con_pl_23], screen_width+0.02)
screening_pl_22.segment((end_screen_finger[con_pl_23][0]+0.04, end_screen_finger[con_pl_23][1]-0.02), con_bar_width)
screening_pl_22.segment(con_pos_dict['east'][-4], con_bar_width)
screening_pl_22.segment((con_pos_dict_2['east'][-4][0]-pad_int_fo_length, con_pos_dict_2['east'][-4][1]), 0.1)
screening_pl_22.segment(con_pos_dict_2['east'][-4], 0.1)

screening_pl_23 = gdstk.FlexPath(screening_gate_dict[con_pl_23].path_spines()[0]
                                  , screen_width+0.005, simple_path=False, layer=layer_b1, ends=(-0.04,0))
screening_pl_23.segment((end_screen_finger[con_pl_23][0]+0.04, end_screen_finger[con_pl_23][1]-0.02), con_bar_width)
screening_pl_23.segment((con_pos_dict['east'][-4][0], con_pos_dict['east'][-4][1]), con_bar_width)

screening_pl_sens_east = gdstk.FlexPath(screening_gate_dict[con_pl_sens_east].path_spines()[0]
                                      , screen_width, simple_path=False, layer=layer_b1)
screening_pl_sens_east.segment((end_screen_finger[con_pl_sens_east][0]+0.01,
                              end_screen_finger[con_pl_sens_east][1]-0.05),
                               con_bar_width)
screening_pl_sens_east.segment((con_pos_dict_2['east'][5][0]-pad_int_fo_length, con_pos_dict_2['east'][5][1]), 0.1)
screening_pl_sens_east.segment(con_pos_dict_2['east'][5], 0.1)


screening_pl_32 = screening_pl_12.copy()
screening_pl_32.mirror((1,0))
screening_pl_32.mirror((0,1))

screening_pl_33 = screening_pl_11.copy()
screening_pl_33.mirror((1,0))
screening_pl_33.mirror((0,1))

screening_pl_21 = screening_pl_22.copy()
screening_pl_21.mirror((1,0))
screening_pl_21.mirror((0,1))

screening_pl_20 = screening_pl_23.copy()
screening_pl_20.mirror((1,0))
screening_pl_20.mirror((0,1))

screening_pl_sens_south = screening_pl_sens_north.copy()
screening_pl_sens_south.mirror((1,0))
screening_pl_sens_south.mirror((0,1))

screening_pl_sens_west = screening_pl_sens_east.copy()
screening_pl_sens_west.mirror((1,0))
screening_pl_sens_west.mirror((0,1))

screening_gate_sensor_cell.add(screening_pl_11,
                               screening_pl_12,
                               screening_pl_20,
                               screening_pl_21,
                               screening_pl_22,
                               screening_pl_23,
                               screening_pl_32,
                               screening_pl_33,
                               screening_pl_sens_north,
                               screening_pl_sens_east,
                               screening_pl_sens_west,
                               screening_pl_sens_south)


# # add all screening gates to cell
screening_gate_cell.add(gdstk.Reference(screening_gate_sensor_cell))

# %%% Main cell with several devices and lithography alignment marks

main = gdstk.Cell("Main")
unit_cell_len = n_row*spacing_hv

main.add(gdstk.cross((0,0), 3*plunger_width, 0.01),
         gdstk.Reference(device, (-unit_cell_len/2*(1-1), -unit_cell_len/2*(1-1)), columns=1, rows=1, spacing = (spacing,spacing)),
         gdstk.Reference(sensor_cell),
         gdstk.Reference(bar_con_cell),
         gdstk.Reference(screening_gate_cell)
         )

main.flatten()


# %%% Save main cell in file
lib = gdstk.Library()
lib.add(main, *main.dependencies(True))
lib.write_gds(os.path.join(path, file_name))
main.write_svg(os.path.join(path, f'343_array_{generation}p{version}gen.svg'))

# %% fanout
file_name = f'343_array_fo_{generation}p{version}gen.gds'

# %%% fanout layers
# CAD layer numbers
pl_layer = 21
bar_top_layer = 5
bar_top_layer_2 = 51
bar_bot_layer = 31
ohm_layer = 3
# 13

fo_north_layers_dict = {pl_layer:[1,6,12],
                        bar_bot_layer:[2,7,10],
                        bar_top_layer:[3,8,13],
                        bar_top_layer_2:[0,5,11],
                        ohm_layer:[4,9]}

fo_east_layers_dict = {pl_layer:[0,6,10,12],
                       bar_top_layer_2:[4,11],
                       bar_top_layer:[1,7],
                       bar_bot_layer:[2,5,9],
                       ohm_layer:[3,8]}

n_north = len(sorted({x for v in fo_north_layers_dict.values() for x in v}))
n_east = len(sorted({x for v in fo_east_layers_dict.values() for x in v}))

fo_south_layers_dict = {pl_layer:n_north-1-np.array(fo_north_layers_dict[pl_layer]),
                        bar_bot_layer:n_north-1-np.array(fo_north_layers_dict[bar_bot_layer]),
                        bar_top_layer:n_north-1-np.array(fo_north_layers_dict[bar_top_layer]),
                        bar_top_layer_2:n_north-1-np.array(fo_north_layers_dict[bar_top_layer_2]),
                        ohm_layer:n_north-1-np.array(fo_north_layers_dict[ohm_layer])}

fo_west_layers_dict = {pl_layer:n_east-1-np.array(fo_east_layers_dict[pl_layer]),
                        bar_bot_layer:n_east-1-np.array(fo_east_layers_dict[bar_bot_layer]),
                        bar_top_layer:n_east-1-np.array(fo_east_layers_dict[bar_top_layer]),
                        bar_top_layer_2:n_east-1-np.array(fo_east_layers_dict[bar_top_layer_2]),
                        ohm_layer:n_east-1-np.array(fo_east_layers_dict[ohm_layer])}

n_south = len(sorted({x for v in fo_south_layers_dict.values() for x in v}))
n_west = len(sorted({x for v in fo_west_layers_dict.values() for x in v}))

# %%% Generate fanout
fog = fanout_generator()

fog.fo_width = 0.1
fog.bondpad_spacing = 40

fog.fo_spacing_north = 2*(x_max2-edge_space2)/(n_north-1)
fog.fo_spacing_south = 2*(x_max2-edge_space2)/(n_south-1)
fog.fo_spacing_east = 2*(y_max2-edge_space2)/(n_east-1)
fog.fo_spacing_west = 2*(y_max2-edge_space2)/(n_west-1)
fog.fo_yoffset_east = 0
fog.fo_yoffset_west = 0
fog.fo2device_buffer = 0

fog.layers_fo_north_dict = fo_north_layers_dict
fog.layers_fo_east_dict = fo_east_layers_dict
fog.layers_fo_south_dict = fo_south_layers_dict
fog.layers_fo_west_dict = fo_west_layers_dict

fo_dict = fog.generate_fo_points(main, xy_max=[x_max2,y_max2])

main_fo = main.copy('main_fo')
fog.generate_fan(main_fo, fo_dict['all'], only_coarse= True, buff = [500,0], ohm_SiN_spacing = 100, ohm_layer=ohm_layer)

# %%% Visual aid
try:
    layout = gdstk.read_rawcells(os.path.join(path, "chip_layout.gds"))
    main_fo.add(gdstk.Reference(layout['TOP']))
    main_fo.flatten()
except:
    pass

# %%% Write and save

lib = gdstk.Library()
lib.add(main_fo, *main_fo.dependencies(True))
lib.write_gds(os.path.join(path, file_name))
main.write_svg(os.path.join(path, f'343_array_fo_{generation}p{version}gen.svg'))
