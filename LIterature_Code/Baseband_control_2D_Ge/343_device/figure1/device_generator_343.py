#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 11:42:00 2020

@author: nicohendrickx, valentinjohn
"""

#%%
import gdstk
import numpy as np

def rot_mat(theta):
    return np.array([[np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]])

def gen_poly(n, sp = None):
    mat = rot_mat(2*np.pi/n)
    if sp is None:
        sp = [np.cos(np.pi/n), np.sin(np.pi/n)]
    poly = [sp]
    for k in range(n):
        poly.append(np.dot(mat, poly[-1]))
    return [tuple(co) for co in poly]

def line_intersection(line1, line2):
    xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
    ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])

    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)
    if div == 0:
       raise Exception('lines do not intersect')

    d = (det(*line1), det(*line2))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    return x, y

def centroid(vertexes):
     _x_list = [vertex [0] for vertex in vertexes]
     _y_list = [vertex [1] for vertex in vertexes]
     _len = len(vertexes)
     _x = sum(_x_list) / _len
     _y = sum(_y_list) / _len
     return(_x, _y)

def connector(c_start, c_stop, y_int, phi1, con_width=0.03, layer=0,
              con_type1_vert=True, con_type2_vert=False):
    con_width_x1 = con_width/np.sin(phi1)
    
    (x0,y0) = c_start
    (x1,y1) = (c_start[0]-con_width_x1, c_start[1])
    if con_type1_vert:
        (x1,y1) = (c_start[0], c_start[1]-con_width_x1)
    (x3,y3) = (x0-np.tan(phi1)*(y_int-(y0)) - con_width_x1, y_int)
    
    phi2 = np.arctan((y_int-y1)/(x3+con_width_x1-x1))
    con_width_x2 = con_width/np.sin(phi2)

    (x4,y4) = c_stop
    (x5,y5) = (c_stop[0]- con_width_x2, c_stop[1])
    if con_type2_vert:
        (x5,y5) = (c_stop[0], c_stop[1]- con_width_x2)
    
    slope1 = (y3-y1)/(x3-x1)
    slope2 = (y4-y3)/(x4-x3)
    (x2,y2) = line_intersection(((x0,y0),(x0+1,y0+slope1)),
                                ((x5,y5), (x5+1,y5+slope2)))
    
    con = gdstk.Polygon([(x0,y0), (x1,y1),
                         (x3,y3), (x4,y4),
                         (x5,y5), (x2,y2)],
                        layer = layer)
    return con

def plunger_cell(spacing = 0.16, plunger_width = 0.1, diag_width = 0.03,
                 layer_b1 = 31, layer_b2 = 5, pl_layer = 21, bar_width = 0.03,
                 rb1_lb2_horizontal=True,
                 con_top=True, con_bottom=True, con_right=True, con_left=True):
    
    spacing_hv = 2**0.5*spacing
    overlap = (plunger_width + 2*bar_width - spacing)/3
    
    pl_cell = gdstk.Cell('plunger_v1') # define cell
    plunger_points = gen_poly(8) # define polygon
    plunger = gdstk.Polygon(plunger_points, layer = pl_layer) # add poly to cell
    pl_cell.add(plunger)
    
    # Diagonal bit of plunger
    diag_seg = gdstk.Polygon([(-diag_width/2, -spacing_hv/2), (diag_width/2, -spacing_hv/2),
                              (diag_width/2, spacing_hv/2), (-diag_width/2, spacing_hv/2)],
                             layer = pl_layer) # define poly
    
    plunger.scale(0.5/np.cos(np.pi/8)*plunger_width) # scale plunger to right size
    pl_cell.add(diag_seg) # add diag poly cell
    
    
    # position of barriers from plunger centre
    x_bar_1 = 2**0.5/4*(plunger_width+bar_width-2*overlap)
    y_bar_1 = 2**0.5/4*(plunger_width+bar_width-2*overlap)
    x_bar_2 = x_bar_1 + 2**0.5/2*(bar_width-overlap)
    y_bar_2 = y_bar_1 + 2**0.5/2*(bar_width-overlap)
    
    # Left top barrier 1
    bar_lt_1 = gdstk.rectangle((-plunger_width/2, -bar_width/2), (plunger_width/2, bar_width/2), layer = layer_b1)
    bar_lt_1.rotate(np.pi/4)
    bar_lt_1.translate(-x_bar_1, y_bar_1)
    pl_cell.add(bar_lt_1)
    
    # Left top barrier 2
    bar_lt_2 = gdstk.rectangle((-plunger_width/2, -bar_width/2), (plunger_width/2, bar_width/2), layer = layer_b2)
    bar_lt_2.rotate(np.pi/4)
    bar_lt_2.translate(-x_bar_2, y_bar_2)
    pl_cell.add(bar_lt_2)
    
    # Right top barrier 1
    bar_rt_1 = gdstk.rectangle((-plunger_width/2, -bar_width/2), (plunger_width/2, bar_width/2), layer = layer_b2)
    bar_rt_1.rotate(-np.pi/4)
    bar_rt_1.translate(x_bar_1, y_bar_1)
    pl_cell.add(bar_rt_1)
    
    # Right top barrier 2
    bar_rt_2 = gdstk.rectangle((-plunger_width/2, -bar_width/2), (plunger_width/2, bar_width/2), layer = layer_b1)
    bar_rt_2.rotate(-np.pi/4)
    bar_rt_2.translate(x_bar_2, y_bar_2)
    pl_cell.add(bar_rt_2)
    
    # Left bottom barrier 1
    bar_lb_1 = gdstk.rectangle((-plunger_width/2, -bar_width/2), (plunger_width/2, bar_width/2), layer = layer_b2)
    bar_lb_1.rotate(-np.pi/4)
    bar_lb_1.translate(-x_bar_1, -y_bar_1)
    pl_cell.add(bar_lb_1)
    
    # Left bottom barrier 2
    bar_lb_2 = gdstk.rectangle((-plunger_width/2, -bar_width/2), (plunger_width/2, bar_width/2), layer = layer_b1)
    bar_lb_2.rotate(-np.pi/4)
    bar_lb_2.translate(-x_bar_2, -y_bar_2)
    pl_cell.add(bar_lb_2)
    
    # Right bottom barrier 1
    bar_rb_1 = gdstk.rectangle((-plunger_width/2, -bar_width/2), (plunger_width/2, bar_width/2), layer = layer_b1)
    bar_rb_1.rotate(np.pi/4)
    bar_rb_1.translate(x_bar_1, -y_bar_1)
    pl_cell.add(bar_rb_1)
    
    # Right bottom barrier 2
    bar_rb_2 = gdstk.rectangle((-plunger_width/2, -bar_width/2), (plunger_width/2, bar_width/2), layer = layer_b2)
    bar_rb_2.rotate(np.pi/4)
    bar_rb_2.translate(x_bar_2, -y_bar_2)
    pl_cell.add(bar_rb_2)
    
    if rb1_lb2_horizontal:
        layer_left2right = layer_b1
        layer_top2bottom = layer_b2
    else:
        layer_left2right = layer_b2
        layer_top2bottom = layer_b1
    
    xy_con = 2**0.5/2*(plunger_width-overlap)
    xy_con_offset = -2**0.5/2*overlap
    
    conn_points = [(0,0), (bar_width,2*bar_width), (0,2*bar_width), (-bar_width, 0)] # define poly points
    conn_points_t = [(0,0), (-bar_width,2*bar_width), (0,2*bar_width), (bar_width, 0)] # define poly points
    
    # connecting the first right bottom barrier with the second left bottom barrier
    conn_point_rb1_lb2 = bar_lb_2.points[1]+2**0.5/2*bar_width*np.array((-1,1))
    conn_points_rb1_lb2 = [bar_rb_1.points[0], bar_rb_1.points[-1], conn_point_rb1_lb2, bar_lb_2.points[1]] # define poly points
    con_rb1_lb2 = gdstk.Polygon(conn_points_rb1_lb2, layer = layer_left2right) # make poly     
    if not rb1_lb2_horizontal:
        con_rb1_lb2.mirror((0,1))
    pl_cell.add(con_rb1_lb2)
    # connecting the first left top barrier with the second right top barrier
    conn_point_lt1_rt2 = bar_rt_2.points[3]+2**0.5/2*bar_width*np.array((1,-1))
    conn_points_lt1_rt2 = [bar_lt_1.points[2], bar_lt_1.points[1], conn_point_lt1_rt2, bar_rt_2.points[3]] # define poly points
    con_lt1_rt2 = gdstk.Polygon(conn_points_lt1_rt2, layer = layer_left2right) # make poly   
    if not rb1_lb2_horizontal:
        con_lt1_rt2.mirror((0,1))
    pl_cell.add(con_lt1_rt2)

    conn_point_lb1_lt2 = bar_lt_2.points[3]+2**0.5/2*bar_width*np.array((1,1))
    conn_points_lb1_lt2 = [bar_lb_1.points[0], bar_lb_1.points[3], conn_point_lb1_lt2, bar_lt_2.points[3]] # define poly points
    con_lb1_lt2 = gdstk.Polygon(conn_points_lb1_lt2, layer = layer_top2bottom) # make poly     
    if not rb1_lb2_horizontal:
        con_lb1_lt2.mirror((1,0))
    pl_cell.add(con_lb1_lt2)
    
    conn_point_rt1_rb2 = bar_rb_2.points[1]+2**0.5/2*bar_width*np.array((-1,-1))
    conn_points_rt1_rb2 = [bar_rt_1.points[2], bar_rt_1.points[1], conn_point_rt1_rb2, bar_rb_2.points[1]] # define poly points
    con_rt1_rb2 = gdstk.Polygon(conn_points_rt1_rb2, layer = layer_top2bottom) # make poly  
    if not rb1_lb2_horizontal:
        con_rt1_rb2.mirror((1,0))
    pl_cell.add(con_rt1_rb2)

    if con_right:
        con_rb1 = con_lb1_lt2.copy()
        con_rb1.layer = layer_left2right
        con_rb1.translate(2**0.5*(plunger_width-overlap),0) # move in right spot
        con_rb1.mirror(con_rb1.points[0], con_rb1.points[1])
        pl_cell.add(con_rb1)
    if con_left:
        con_lt1 = con_rt1_rb2.copy()
        con_lt1.layer = layer_left2right
        con_lt1.translate(-2**0.5*(plunger_width-overlap),0) # move in right spot
        con_lt1.mirror(con_lt1.points[0], con_lt1.points[1])
        pl_cell.add(con_lt1)
    if con_top:
        con_rt1 = con_rb1_lb2.copy()
        con_rt1.layer = layer_top2bottom
        con_rt1.translate(0,2**0.5*(plunger_width-overlap)) # move in right spot
        con_rt1.mirror(con_rt1.points[0], con_rt1.points[1])
        pl_cell.add(con_rt1)
    if con_bottom:
        con_lb1 = con_lt1_rt2.copy()
        con_lb1.layer = layer_top2bottom
        con_lb1.translate(0,-2**0.5*(plunger_width-overlap)) # move in right spot
        con_lb1.mirror(con_lb1.points[0], con_lb1.points[1])
        pl_cell.add(con_lb1)

    return pl_cell

def gen_sensor(orientation_sensor='N', orientation_source='W', orientation_drain='E', 
               pl_layer=21, ohm_layer=3, layer_bar_sen=(31, 31),
               sensor_width=0.16, sensor_sd_width = 0.030,
               sensor_position=(0,0), bar_width=0.04, bar_length=0.05,
               feature_gap=0, plunger_gap_short = 0.035,
               source_position_offset=(0,0),  drain_position_offset=(0,0),
               bar_sou_position_offset=(0,0), bar_dra_position_offset=(0,0),
               gen_sd = True,
               bar_sharp_source = 0, bar_sharp_drain = 0):
    sensor_cell = gdstk.Cell('sensor_'+str(orientation_sensor)) 
    sensor_points = gen_poly(8)
    sensor_pl = gdstk.Polygon(sensor_points, layer = pl_layer)
    sensor_pl.scale(0.5/np.cos(np.pi/8)*sensor_width) # scale plunger to right size
    
    orientation_dict = {'N':(0,1), 'E':(1,0), 'S':(0,-1), 'W':(-1,0),
                        'NE':(2**0.5/2,2**0.5/2), 
                        'SE':(2**0.5/2,-2**0.5/2), 
                        'SW':(-2**0.5/2,-2**0.5/2), 
                        'NW':(-2**0.5/2,2**0.5/2)}
    
    bar_angle_dict = {'N':0, 'E':-np.pi/2, 'S':np.pi, 'W':+np.pi/2,
                      'NE':-np.pi/4, 
                      'SE':np.pi/4, 
                      'SW':3*np.pi/4, 
                      'NW':-3*np.pi/4}
    
    (i,j) = orientation_dict[orientation_source]
    (m,n) = orientation_dict[orientation_drain]
    
    sd_position = ((i*(sensor_width/2+bar_width+sensor_sd_width/2)+source_position_offset[0],
                    j*(sensor_width/2+bar_width+sensor_sd_width/2)+source_position_offset[1]),
                   (m*(sensor_width/2+bar_width+sensor_sd_width/2)+drain_position_offset[0],
                    n*(sensor_width/2+bar_width+sensor_sd_width/2)+drain_position_offset[1]))
    
    bar_position = ((i*(sensor_width/2+bar_width/2-feature_gap)+bar_sou_position_offset[0],
                     j*(sensor_width/2+bar_width/2-feature_gap)+bar_sou_position_offset[1]),
                    (m*(sensor_width/2+bar_width/2-feature_gap)+bar_dra_position_offset[0],
                     n*(sensor_width/2+bar_width/2-feature_gap)+bar_dra_position_offset[1]))
    
    sensor_source = gdstk.Polygon(sensor_points, layer = ohm_layer)
    sensor_source.scale(0.5/np.cos(3.14159265359/8)*sensor_sd_width) # scale plunger to right size
    sensor_source.translate(sd_position[0])
    sensor_drain = gdstk.Polygon(sensor_points, layer = ohm_layer)
    sensor_drain.scale(0.5/np.cos(3.14159265359/8)*sensor_sd_width) # scale plunger to right size
    sensor_drain.translate(sd_position[1])
     
    bar_sen_source = gdstk.Polygon([(-bar_width/2, -bar_length/2),
                                    (-bar_width/2, bar_length/2),
                                    (bar_width/2, bar_length/2),
                                    (bar_width/2, -bar_length/2-bar_sharp_source)
                                    ],
                                        layer = layer_bar_sen[0])
    
    bar_sen_source.rotate(bar_angle_dict[orientation_sensor])
    bar_sen_source.translate(bar_position[0][0], bar_position[0][1])
    
    
    bar_sen_drain = gdstk.Polygon([(-bar_width/2, -bar_length/2-bar_sharp_drain),
                                    (-bar_width/2, bar_length/2),
                                    (bar_width/2, bar_length/2),
                                    (bar_width/2, -bar_length/2)
                                    ],
                                        layer = layer_bar_sen[1])
    bar_sen_drain.rotate(bar_angle_dict[orientation_sensor])
    bar_sen_drain.translate(bar_position[1][0], bar_position[1][1])
   
    sensor_cell.add(sensor_pl)
    
    if gen_sd:
        sensor_cell.add(
               sensor_source,
               sensor_drain,
               )
    sensor_cell.add(bar_sen_source,
               bar_sen_drain
               )
    

    # if bar_sharp_drain:
    #     bar_sen_drain_sharp = gdstk.Polygon([(bar_width/2, 0),
    #                                           (-bar_width/2, 0),
    #                                           (-bar_width/2+bar_sou_position_offset[0], -bar_sou_position_offset[1])],
    #                                           layer = layer_bar_sen)
    #     bar_sen_drain_sharp.rotate(bar_angle_dict[orientation_drain])
    #     bar_sen_drain_sharp.translate(bar_position[1][0], bar_position[1][1])
    #     bar_sen_drain_sharp.translate(-bar_length/2*orientation_dict[orientation_sensor][0], -bar_length/2*orientation_dict[orientation_sensor][1])
    #     sensor_cell.add(bar_sen_drain_sharp)
        
    return sensor_cell, sd_position

def flatten(t):
    return [item for sublist in t for item in sublist]

def fo_layers(fo_layers_dict):
    fo_layer = []
    for key in fo_layers_dict.keys():
        fo_layer.append(fo_layers_dict[key])
    flatten(fo_layer)
    fo_layer = flatten(fo_layer)
    fo_layer.sort()
    n = len(fo_layer)
    if not fo_layer == list(range(n)):
        print('Assignemnts of layers for fanout lines may be incomplete or ambiguous')
    # check if flatten(fo_pl_layer, fo_bar_top_layer, fo_bar_bot_layer) is unique
    
    fo_layer_inverted = np.zeros(n, dtype=int)
    for key in fo_layers_dict.keys():
        for i in fo_layers_dict[key]:
            fo_layer_inverted[i] = key
    
    return list(map(int, fo_layer_inverted))

from dataclasses import dataclass
@dataclass
class fo_point():
    p1: tuple
    p2: tuple
    layer: int
    side: str
    index: int = None

    def __post_init__(self):
        self.sort_points()
    
    def sort_points(self):
        x1 = self.p1[0]
        x2 = self.p2[0]
        y1 = self.p1[1]
        y2 = self.p2[1]
        if self.side == 'north' or self.side == 'south':
            if x2 < x1:
                self.p1 = (x2, y2)
                self.p2 = (x1, y1)
        elif self.side == 'west' or self.side == 'east':
            if y2 < y1:
                self.p1 = (x2, y2)
                self.p2 = (x1, y1)

class fanout_generator():
    def __init__(self):
        self.layer_dict = None
        self.fo_width = 0.03
        self.fo_spacing_north = 0.07
        self.fo_spacing_east = 0.07
        self.fo_spacing_south = 0.07
        self.fo_spacing_west = 0.07
        self.fo_xoffset_north = 0
        self.fo_xoffset_south = 0
        self.fo_yoffset_east = 0
        self.fo_yoffset_west = 0
        self.fo2device_buffer = 0.2
        self.layers_fo_north_dict = {0:list(np.arange(16,dtype=int))}
        self.layers_fo_east_dict = {0:list(np.arange(16,dtype=int))}
        self.layers_fo_south_dict = {0:list(np.arange(16,dtype=int))}
        self.layers_fo_west_dict = {0:list(np.arange(16,dtype=int))}
        self.chip_width = 4000
        self.margin = 300        
        self.bondpad_spacing = 100
        self.bondpad_length = 400        
        self.fan_width = 25
        self.int_spacing = 1
        self.int_width = 2
        self.use_width = self.chip_width - 2 * self.margin - 2 * self.bondpad_length - self.bondpad_spacing/2
        self.use_height = self.chip_width - 2 * self.margin   

    def course_layer(self, layer):
        if self.layer_dict:
            return self.layer_dict[layer]
        else:
            return layer + 1
    
    def generate_fo_points(self, device, xy_max=None):
        fo_points = list()
        fo_points_north = list()
        fo_points_south = list()
        fo_points_west = list()
        fo_points_east = list()
        
        buffer = self.fo2device_buffer
        # determine max and min values of device automatically
        if xy_max == None:
            x_max = 0
            y_max = 0
            x_min = 0
            y_min = 0
            for polygon in device.polygons:
                for point in polygon.points:
                    x = point[0]
                    y = point[1]
                    if x > x_max:
                        x_max = x
                    if y > y_max:
                        y_max = y
                    if x < x_min:
                        x_min = x
                    if y < y_min:
                        y_min = y
        
        else:
            x_max = xy_max[0]
            y_max = xy_max[1]
            x_min = -x_max
            y_min = -y_max
        print(x_max, y_max)

        layers_fo_north = fo_layers(self.layers_fo_north_dict)
        layers_fo_east = fo_layers(self.layers_fo_east_dict)
        layers_fo_south = fo_layers(self.layers_fo_south_dict)
        layers_fo_west = fo_layers(self.layers_fo_west_dict)

        n_north = len(layers_fo_north)
        n_east = len(layers_fo_east)
        n_south = len(layers_fo_south)
        n_west = len(layers_fo_west)

        for n in range(n_north):
            h_y1 = y_max + buffer
            h_y2 = y_max + buffer
            h_x1 = -((n_north-1)/2-n)*self.fo_spacing_north - self.fo_width/2 + self.fo_xoffset_north
            h_x2 = -((n_north-1)/2-n)*self.fo_spacing_north + self.fo_width/2 + self.fo_xoffset_north
            fo_hb = fo_point((h_x1, h_y1), (h_x2, h_y2),
                             layers_fo_north[n], 'north')
            fo_points_north.append(fo_hb)   
            fo_points.append(fo_hb)

        for n in range(n_east):
            h_y1 = -((n_east-1)/2-n)*self.fo_spacing_east - self.fo_width/2 + self.fo_yoffset_east
            h_y2 = -((n_east-1)/2-n)*self.fo_spacing_east + self.fo_width/2 + self.fo_yoffset_east
            h_x1 = x_max + buffer
            h_x2 = x_max + buffer
            fo_hb = fo_point((h_x1, h_y1), (h_x2, h_y2),
                             layers_fo_east[n], 'east')
            fo_points_east.append(fo_hb)   
            fo_points.append(fo_hb)

        for n in range(n_south):
            h_y1 = y_min - buffer
            h_y2 = y_min - buffer
            h_x1 = -((n_south-1)/2-n)*self.fo_spacing_south - self.fo_width/2 + self.fo_xoffset_south
            h_x2 = -((n_south-1)/2-n)*self.fo_spacing_south + self.fo_width/2 + self.fo_xoffset_south
            fo_hb = fo_point((h_x1, h_y1), (h_x2, h_y2),
                             layers_fo_south[n], 'south')
            fo_points_south.append(fo_hb)   
            fo_points.append(fo_hb)

        for n in range(n_west):
            h_y1 = -((n_west-1)/2-n)*self.fo_spacing_west - self.fo_width/2 + self.fo_yoffset_west
            h_y2 = -((n_west-1)/2-n)*self.fo_spacing_west + self.fo_width/2 + self.fo_yoffset_west
            h_x1 = x_min - buffer
            h_x2 = x_min - buffer
            fo_hb = fo_point((h_x1, h_y1), (h_x2, h_y2),
                             layers_fo_west[n], 'west')
            fo_points_west.append(fo_hb)   
            fo_points.append(fo_hb) 

        return {'all':fo_points,
                'north':fo_points_north,
                'east':fo_points_east,
                'south':fo_points_south,
                'west':fo_points_west}

    def generate_fan(self, device, fo_points, only_coarse = False, buff = None, ohm_layer = 3, ohm_SiN_spacing = 50):
        if buff is None:
            buff = [0,0]
        sides = ['west', 'north', 'east', 'south']
        for side in sides:
            # define sorting axis
            if side == 'west' or side == 'east':
                sort_ind = 1
            elif side == 'north' or side == 'south':
                sort_ind = 0
            
            # grab relevant points
            relevant_fo = [fo for fo in fo_points if fo.side == side]
            
            # sort to avoid crossing fanouts
            p1_inds = [rfo.p1[sort_ind] for rfo in relevant_fo]
            indices = np.array(p1_inds).argsort().argsort()
            print(f'{side}: {[rfo.p1[sort_ind] for rfo in relevant_fo]}')
            print(f'{side}: {indices}')
            #calculate int_sizes
            tot_int_size = 2*len(relevant_fo) * self.int_spacing
            bondpad_length = self.bondpad_length
            # set index to fo object and generate fanout
            
            j = - 1
            for (fo, ind) in zip(relevant_fo, indices):
                fo.index = ind
                k = fo.index
                ns = 1 if side == 'north' or side == 'east' else -1
                int_points = [(-tot_int_size/2 + 2*(k+0.25)*self.int_spacing, ns*tot_int_size/2),
                            (-tot_int_size/2 + 2*(k+0.25)*self.int_spacing, ns*tot_int_size/2 + ns*self.int_width),
                            (-tot_int_size/2 + 2*(k+0.75)*self.int_spacing, ns*tot_int_size/2 + ns*self.int_width),
                            (-tot_int_size/2 + 2*(k+0.75)*self.int_spacing, ns*tot_int_size/2)]
                
                if side == 'east' or side == 'west':
                    int_points = [ip[::-1] for ip in int_points]
                
                if only_coarse:
                    int_points = [fo.p2, fo.p1, *int_points]
                    int_piece = gdstk.Polygon(int_points, layer = fo.layer)
                    device.add(int_piece)
                    
                    if side == 'east' or side == 'west':
                        uw = self.use_width - buff[0]
                    else:    
                        uw = self.use_width - buff[1]           
                    uh = self.use_height
                    bp_width = (uw - (len(relevant_fo) - 1) * self.bondpad_spacing)/len(relevant_fo)
                    bp_period = bp_width + self.bondpad_spacing
                    
                    x1 = (-uw/2 + k * bp_period + (bp_width - self.fan_width)/2)/4
                    y1 = (ns * uh/2 - ns*bondpad_length)/2
                    
                    x2 = -uw/2 + k * bp_period + (bp_width - self.fan_width)/2
                    y2 = ns * uh/2 - ns*bondpad_length
                    
                    x3 = -uw/2 + k * bp_period
                    y3 = ns * uh/2 - ns*bondpad_length
                    
                    x4 = -uw/2 + k * bp_period
                    y4 = ns * uh/2
                    
                    x5 = -uw/2 + k * bp_period + bp_width
                    y5 = ns * uh/2
                    
                    x6 = -uw/2 + k * bp_period + bp_width
                    y6 = ns * uh/2 - ns*bondpad_length
                    
                    x7 = -uw/2 + k * bp_period + (bp_width + self.fan_width)/2
                    y7 = ns * uh/2 - ns*bondpad_length
                    
                    x8 = (-uw/2 + k * bp_period + (bp_width + self.fan_width)/2)/4
                    y8 = y1
                    
                    if fo.layer != ohm_layer:
                        fan_points = [(x1, y1),
                                      (x2, y2),
                                      (x3, y3),
                                      (x4, y4),                              
                                      (x5, y5),                          
                                      (x6, y6),
                                      (x7, y7),
                                      (x8, y8)]
                        print(fan_points)
                    else: # this needs to be changed for ohmics !!!!!!!!!!!!!!!
                        xs2 = -uw/2 + k * bp_period + bp_width/2
                        xs1 = xs2/4
                        ys2 = ns * uh/2 - ns*bondpad_length
                        ys1 = ys2/2
                        slope = (xs2-xs1)/(ys2-ys1)
                        delta_y = bondpad_length
                        delta_x = delta_y/slope
                        # print('delta_x = '+str(delta_x))
                        # print('x1= '+str(x1))
                        # print('x2 = '+str(x2))
                        bp_ohmic_width_end = bp_width+2*self.bondpad_spacing
                        bp_ohmic_width = bp_width*(bondpad_length/(ys2-ys1))
                        
                        c3 = (x3- ((ns*(1*ohm_SiN_spacing+bondpad_length)*slope)) - (bp_ohmic_width-bp_width)/2,
                              y3-ns*(bondpad_length+ohm_SiN_spacing))
                        c6 = (x6- ((ns*(1*ohm_SiN_spacing+bondpad_length)*slope)) + (bp_ohmic_width-bp_width)/2,
                              y6- ns*(bondpad_length+ohm_SiN_spacing))
                        
                        if ns != 1:
                            (c3,c6) = (c6,c3)
                        
                        fan_points = [(x1,y1),
                                      (x2- ns*(ohm_SiN_spacing+bondpad_length)*slope, y2-ns*(bondpad_length+ohm_SiN_spacing)),
                                      c3,
                                      (x4 - (ns*ohm_SiN_spacing*slope + (bp_ohmic_width_end-bp_width)/2), y4-ns*(bondpad_length+ohm_SiN_spacing)),
                                      (x5- (ns*ohm_SiN_spacing*slope - (bp_ohmic_width_end-bp_width)/2), y5-ns*(bondpad_length+ohm_SiN_spacing)),
                                      c6,
                                      (x7- ns*(ohm_SiN_spacing+bondpad_length)*slope, y7-ns*(bondpad_length+ohm_SiN_spacing)),
                                      (x8,y8)
                                      ]

                    
                    if side == 'east' or side == 'west':
                        fan_points = [fp[::-1] for fp in fan_points]
                    
                    fan_points = [*int_points[4:6], *int_points[2:4], *fan_points]
                    fan_piece = gdstk.Polygon(fan_points, layer = self.course_layer(fo.layer))
                    device.add(fan_piece)
                else:                    
                    bp_width = (self.use_width - (len(relevant_fo) - 1) * self.bondpad_spacing)/len(relevant_fo)
                    bp_period = bp_width + self.bondpad_spacing
                    fan_points = [((-self.use_width/2 + k * bp_period + (bp_width - self.fan_width)/2)/4, (ns * self.use_height/2 - ns*bondpad_length)/2),
                                  (-self.use_width/2 + k * bp_period + (bp_width - self.fan_width)/2, ns * self.use_height/2 - ns*bondpad_length),
                                  (-self.use_width/2 + k * bp_period, ns * self.use_height/2 - ns*bondpad_length),
                                  (-self.use_width/2 + k * bp_period, ns * self.use_height/2),                              
                                  (-self.use_width/2 + k * bp_period + bp_width, ns * self.use_height/2),                          
                                  (-self.use_width/2 + k * bp_period + bp_width, ns * self.use_height/2 - ns*bondpad_length),
                                  (-self.use_width/2 + k * bp_period + (bp_width + self.fan_width)/2, ns * self.use_height/2 - ns*bondpad_length),
                                  ((-self.use_width/2 + k * bp_period + (bp_width + self.fan_width)/2)/4, (ns * self.use_height/2 - ns*bondpad_length)/2)]
                    
                    if side == 'east' or side == 'west':
                        fan_points = [fp[::-1] for fp in fan_points]
                    
                    fan_points = [fo.p2, fo.p1, *fan_points]
                    fan_piece = gdstk.Polygon(fan_points, layer = self.course_layer(fo.layer))
                    device.add(fan_piece)
                
                    