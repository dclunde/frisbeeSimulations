# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 11:02:52 2016

@author: dlunde
"""

import math
import random
import numpy as np
import pylab
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider, Button,RadioButtons
from scipy.linalg import norm
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
import matplotlib as plotlib
#plotlib.use('TkAgg') #Qt4Agg
import sys
import matplotlib.patches as patches
import operator
from matplotlib.animation import FuncAnimation
#from subprocess import call
import os
import time

if (sys.version_info > (3, 0)):
     # Python 3 code in this block
     import matplotlib.pyplot as plt,mpld3
     from bokeh import mpl
     from bokeh.plotting import save
     #import bokeh     
     print "Boken uploaded"

else:
     # Python 2 code in this block
    import matplotlib.pyplot as plt


#==============================================================================
#global x,y,z,vx,vy,vz,deltavx,deltavy,deltavz,delta_t,mass,\
#    delta_omega_x,delta_omega_y,delta_omega_z,Ixx,Iyy,Izz,omega_x,omega_y,omega_z,\
#    alpha,phi,mat_x,mat_y,mat_z,mat_vz,mat_vy,mat_vz,mat_a,mat_p,mat_t,mat_wx,mat_wy,mat_wz
    
# Flight Variables
alpha    = 8        # Angle of attack
alpha_0  = -4       # Initial Angle of Attack according to frisbee
phi      = 0        # Z rotation Angle
wind     = [0,0,0]  # 3D vector of wind

# Initial Variables
x     = 0.0   # X initial position
y     = 1.0   # Y initial position
z     = 0.0   # Z initial position
vx    = 14    # velocity * math.cos(math.radians(alpha-alpha_0))   # Initial velocity in X position
vy    = 3     # velocity * math.sin(math.radians(alpha-alpha_0))   # Initial velocity in Y position
vz    = 0     # Initial Z velocity

#x_0     = 0.0   # X initial position
#y_0     = 1.0   # Y initial position
#z_0     = 0.0   # Z initial position
#vx_0    = 14    # velocity * math.cos(math.radians(alpha-alpha_0))   # Initial velocity in X position
#vy_0    = 3     # velocity * math.sin(math.radians(alpha-alpha_0))   # Initial velocity in Y position
#vz_0    = 0     # Initial Z velocity

omega_x   = 0        # Initial speed of rotation in x direction (rot/s)
omega_y   = 10      # Initial speed of rotation in y direction (rot/s)
omega_z   = 0        # Initial speed of rotation in z direction (rot/s)

#==============================================================================
#Vector
arrows = [0]*4
frisbee_parts =  [0]*3

delta_t = 0.01  # Change in time

# Frisbee Variables
mass = 0.175                      # Mass of Frisbee in kg in the AUDL
diameter = 0.27305
area = math.pi*pow(diameter/2.0,2) # Area (m^2) of normal frisbee in the AUDL

# Working Variables
#x    = 0.0    # X position
#y    = 0.0    # Y position
#z    = 0.0    # Z position
#vx   = 0.0    # Velocity in X position
#vy   = 0.0    # Velocity in Y position
#vz   = 0.0    # Velocity in Z position
cl   = 0.0    # Cooeficient of Life
cd   = 0.0    # Coeficient of Drag
cl_w = 0.0    # Wind Lift
cd_w = 0.0    # Wind Drag
deltavx = 0.0 # Change in X Velocity
deltavy = 0.0 # Change in Y Velocity
deltavz = 0.0 # Change in Z Velocity
delta_omega_x   = 0 # Initial acceleration of rotation in x direction
delta_omega_y   = 0 # Initial acceleration of rotation in y direction
delta_omega_z   = 0 # Initial acceleration of rotation in z direction
g    = -9.81  # Gravity constant

#Matrix Creation
total_t = 20                         # total_t/delta_t makes size total_t seconds
mat_x  = [0.0]*int(total_t/delta_t)  # initiate X Distance Matrix
mat_y  = [0.0]*int(total_t/delta_t)  # initiate Y Distance Matrix
mat_z  = [0.0]*int(total_t/delta_t)  # initiate Z Distance Matrix

mat_vx = [0.0]*int(total_t/delta_t)  # initiate X Velocity Matrix
mat_vy = [0.0]*int(total_t/delta_t)  # initiate Y Velocity Matrix
mat_vz = [0.0]*int(total_t/delta_t)  # initiate Z Velocity Matrix

mat_a  = [0.0]*int(total_t/delta_t)  # initiate Angle of Attack Matrix
mat_p  = [0.0]*int(total_t/delta_t)  # initiate Phi Matrix

mat_t  = [0.0]*int(total_t/delta_t)  # initiate Time Matrix

mat_wx  = [0.0]*int(total_t/delta_t) # initiate Omega in X direction Matrix
mat_wy  = [0.0]*int(total_t/delta_t) # initiate Omega in Y direction  Matrix
mat_wz  = [0.0]*int(total_t/delta_t) # initiate Omega in Z direction Matrix

mat_KE   = [0.0]*int(total_t/delta_t) # initiate Kinetic Energy Total
mat_PE   = [0.0]*int(total_t/delta_t) # initiate Potential Energy Total
mat_AE   = [0.0]*int(total_t/delta_t) # initiate Angular Kinetic Energy Total

mat_vxE = [0.0]*int(total_t/delta_t) # initiate Energy in X direction Matrix
mat_vyE = [0.0]*int(total_t/delta_t) # initiate Energy in X direction Matrix
mat_vzE = [0.0]*int(total_t/delta_t) # initiate Energy in X direction Matrix

mat_wxE = [0.0]*int(total_t/delta_t) # initiate Angular Energy in X direction Matrix
mat_wyE = [0.0]*int(total_t/delta_t) # initiate Angular Energy in X direction Matrix
mat_wzE = [0.0]*int(total_t/delta_t) # initiate Angular Energy in X direction Matrix

#Create Vectors
velocity_vec         = [[0,0,0]]*int(total_t/delta_t)
frisbee_vec          = [[0,0,0]]*int(total_t/delta_t)
velocity_on_frisbee  = [[0,0,0]]*int(total_t/delta_t)
lift_vec             = [[0,0,0]]*int(total_t/delta_t)
wind_on_frisbee      = [[0,0,0]]*int(total_t/delta_t)
lift_wind_vec        = [[0,0,0]]*int(total_t/delta_t)

def grab_from_file():
#    global x_0,y_0,z_0,vx_0,vy_0,vz_0,omega_x,omega_y,omega_z,alpha,phi
    global x,y,z,vx,vy,vz,omega_x,omega_y,omega_z,alpha,phi
    ff = open('frisbee_input.txt','r')
    input_line = ff.readline()
    input_line = ff.readline()
    ff.close()
    input_line = input_line.strip('\n')
    input_mat = [float(fff) for fff in input_line.split(',')]
    if len(input_mat) != 11:
        print 'Wrong input file'
        return
#    x_0,y_0,z_0,vx_0,vy_0,vz_0,omega_x,omega_y,omega_z,alpha,phi = input_mat
    x,y,z,vx,vy,vz,omega_x,omega_y,omega_z,alpha,phi = input_mat
    print input_mat,phi,"<<PHI!"
    return

def edit_initial_conditions():
    global x,y,z,vx,vy,vz,omega_x,omega_y,omega_z,alpha,phi,alpha_0
    ff = open('frisbee_input.txt','r')
    input_line1 = ff.readline()
    input_line = ff.readline()
    ff.close()
    input_line = input_line.strip('\n')        
    input_line1 = input_line1.strip('\n')        
    input_mat = [float(item) for item in input_line.split(',')]
    labels = [item.strip(" ") for item in input_line1.split(',')]
    if len(input_mat) != 11:
        print 'Wrong input file'
        return
    
    editing = True
    while editing:    
        for i in range(0,np.size(input_mat),1):
            print i,"-",labels[i],"-",input_mat[i]
        test = 0
        test = input("Enter number on left to edit value\rEnter -2 to reset numbers\rEnter -1 to run Simulation :")
        
        if test in [0,1,2,3,4,5,6,7,8,9,10]:
            string = "Enter new value for "+labels[test]+" : "
            new_val = input(string)
            input_mat[test] = int(new_val)
            print new_val
        elif test == -2:
            input_mat = [0,1,0,14,3,0,0,10,0,8,0]
        elif test == -1:
            editing = False
        else:
            input_mat = [0,1,0,14,3,0,0,10,0,8,0]
            print "Not valid nubmer, using a normal backhand"
            editing = False
        
    x,y,z,vx,vy,vz,omega_x,omega_y,omega_z,alpha,phi = input_mat
    
    fff = open('frisbee_input.txt','w')
    for item in range(0,10,1):
        print >> fff,labels[item],",",

    print >> fff,labels[10],'\n',
    
    for item in range(0,10,1):
        print >> fff,input_mat[item],",",

    print >> fff,input_mat[10]
    fff.close()  
    return

def choose_throw_type():
    global x,y,z,vx,vy,vz,omega_x,omega_y,omega_z,alpha,phi,alpha_0
    test = 0
    test = input("0 - Forehand \r1 - Backhand \r2 - Bounce \r3 - High Drop \r4 - Backhand left \rChoose throw type by typing number on left :")

    if test == 0:
        # Forehand
        alpha    = 8        # Angle of attack
        alpha_0  = -4       # Initial Angle of Attack according to frisbee
        phi      = 0        # Z rotation Angle
        
        # Initial Variables
        x     = 0.0   # X initial position
        y     = 1.0   # Y initial position
        z     = 0.0   # Z initial position
        vx    = 14    # velocity * math.cos(math.radians(alpha-alpha_0))   # Initial velocity in X position
        vy    = 3     # velocity * math.sin(math.radians(alpha-alpha_0))   # Initial velocity in Y position
        vz    = 0     # Initial Z velocity
        
        omega_x   = 0        # Initial speed of rotation in x direction (rot/s)
        omega_y   = 10      # Initial speed of rotation in y direction (rot/s)
        omega_z   = 0        # Initial speed of rotation in z direction (rot/s)
    
    elif test == 1:
        # Backhand
        alpha    = 8        # Angle of attack
        alpha_0  = -4       # Initial Angle of Attack according to frisbee
        phi      = 0        # Z rotation Angle
        
        # Initial Variables
        x     = 0.0   # X initial position
        y     = 1.0   # Y initial position
        z     = 0.0   # Z initial position
        vx    = 14    # velocity * math.cos(math.radians(alpha-alpha_0))   # Initial velocity in X position
        vy    = 3     # velocity * math.sin(math.radians(alpha-alpha_0))   # Initial velocity in Y position
        vz    = 0     # Initial Z velocity
        
        omega_x   = 0        # Initial speed of rotation in x direction (rot/s)
        omega_y   = -10      # Initial speed of rotation in y direction (rot/s)
        omega_z   = 0        # Initial speed of rotation in z direction (rot/s)
        
    elif test == 2:
        # Bounce
        alpha    = -8        # Angle of attack
        alpha_0  = -10       # Initial Angle of Attack according to frisbee
        phi      = 0        # Z rotation Angle
        
        # Initial Variables
        x     = 0.0   # X initial position
        y     = 1.0   # Y initial position
        z     = 0.0   # Z initial position
        vx    = 10    # velocity * math.cos(math.radians(alpha-alpha_0))   # Initial velocity in X position
        vy    = -7     # velocity * math.sin(math.radians(alpha-alpha_0))   # Initial velocity in Y position
        vz    = 0     # Initial Z velocity
        
        omega_x   = 0        # Initial speed of rotation in x direction (rot/s)
        omega_y   = 10      # Initial speed of rotation in y direction (rot/s)
        omega_z   = 0        # Initial speed of rotation in z direction (rot/s)
#    elif test == 3:
#        
        
    elif test == 3:
        # High drop        
        alpha    = -4        # Angle of attack
        alpha_0  = -4       # Initial Angle of Attack according to frisbee
        phi      = 0        # Z rotation Angle
        
        # Initial Variables
        x     = 0.0   # X initial position
        y     = 5.0   # Y initial position
        z     = 0.0   # Z initial position
        vx    = 1    # velocity * math.cos(math.radians(alpha-alpha_0))   # Initial velocity in X position
        vy    = 0     # velocity * math.sin(math.radians(alpha-alpha_0))   # Initial velocity in Y position
        vz    = 0     # Initial Z velocity
        
        omega_x   = 0        # Initial speed of rotation in x direction (rot/s)
        omega_y   = 10      # Initial speed of rotation in y direction (rot/s)
        omega_z   = 0        # Initial speed of rotation in z direction (rot/s)
        
    elif test ==4:
        # Backhand left
        alpha    = 8        # Angle of attack
        alpha_0  = -4       # Initial Angle of Attack according to frisbee
        phi      = 10        # Z rotation Angle
        
        # Initial Variables
        x     = 0.0   # X initial position
        y     = 1.0   # Y initial position
        z     = 0.0   # Z initial position
        vx    = 14    # velocity * math.cos(math.radians(alpha-alpha_0))   # Initial velocity in X position
        vy    = 3     # velocity * math.sin(math.radians(alpha-alpha_0))   # Initial velocity in Y position
        vz    = 0     # Initial Z velocity
        
        omega_x   = 0        # Initial speed of rotation in x direction (rot/s)
        omega_y   = 10      # Initial speed of rotation in y direction (rot/s)
        omega_z   = 0        # Initial speed of rotation in z direction (rot/s)
    else:
        print "Not a possible throw, using forehand"
        # Forehand
        alpha    = 8        # Angle of attack
        alpha_0  = -4       # Initial Angle of Attack according to frisbee
        phi      = 0        # Z rotation Angle
        
        # Initial Variables
        x     = 0.0   # X initial position
        y     = 1.0   # Y initial position
        z     = 0.0   # Z initial position
        vx    = 14    # velocity * math.cos(math.radians(alpha-alpha_0))   # Initial velocity in X position
        vy    = 3     # velocity * math.sin(math.radians(alpha-alpha_0))   # Initial velocity in Y position
        vz    = 0     # Initial Z velocity
        
        omega_x   = 0        # Initial speed of rotation in x direction (rot/s)
        omega_y   = 10      # Initial speed of rotation in y direction (rot/s)
        omega_z   = 0        # Initial speed of rotation in z direction (rot/s)
        
    return
    

#CONSTANTS
#if alpha == 0:
CLO  = 0.15    # Lift coefficient at alpha = 0
CDO  = 0.08   # The drag coefficent at alpha = 0
rho  = 1.225      # Density of air in kg/m^3
CLA  = 1.4        # Lift coefficient dependent on alpha
CDA  = 2.72       # The drag coefficent dependent on alpha
B_mag = 0.0001        # Magnus force cooeficient
W_mag = -0.062      # amount of anglular velocity lost per time
CRr = 0.014
CRp = -0.0055
CM0 = -0.08
CMa = 0.43
CMq = -0.005
CNr = -0.0000071
Iyy = 0.002352
Ixx = 0.001219
Izz = Ixx

time_slider = 0
playing = False

def which_planet():
    global rho,g
    test = 0
    test = input("0 - Earth \r1 - Moon \r2 - Mars\r Choose planet :\n")

    if test == 0:
        # Earth
        rho  = 1.225
        g    = -9.81
    elif test == 1:
        # Moon
        rho  = 1.225*10**(-15)
        g    = -1.622
    elif test == 2:
        # Mars
        rho  = 1.225 * 0.006
        g    = -3.711
    else:
        print "Not a possible planet, using Earth"
        test = 0
    return
        
def initialize_variables(overwrite=[]):
#    global mat_x,mat_y,mat_z,mat_vx,mat_vy,mat_vz,mat_a,mat_wx,mat_wy,mat_wz, \
#        x,y,z,vx,vy,vz,x_0,y_0,z_0,vx_0,vy_0,vz_0,alpha,phi,omega_x,omega_y,omega_z 
    global mat_x,mat_y,mat_z,mat_vx,mat_vy,mat_vz,mat_a,mat_p,mat_wx,mat_wy,mat_wz, \
        x,y,z,vx,vy,vz,omega_x,omega_y,omega_z,alpha,phi
    
#    phi = phi + 90
#    omega_x = omega_x * 2 * math.pi
#    omega_y = omega_y * 2 * math.pi
#    omega_z = omega_z * 2 * math.pi
    
#    mat_x[0]  = x_0
#    mat_y[0]  = y_0
#    mat_z[0]  = z_0
#    mat_vx[0] = vx_0
#    mat_vy[0] = vy_0
#    mat_vz[0] = vz_0
    if not overwrite==[]:
        x,y,z,vx,vy,vz,omega_x,omega_y,omega_z,alpha,phi=overwrite
    

    mat_x[0]  = x
    mat_y[0]  = y
    mat_z[0]  = z
    mat_vx[0] = vx
    mat_vy[0] = vy
    mat_vz[0] = vz    
    mat_a[0]  = alpha
    mat_p[0]  = phi + 90
    mat_wx[0] = omega_x * 2 * math.pi
    mat_wy[0] = omega_y * 2 * math.pi
    mat_wz[0] = omega_z * 2 * math.pi
    

#    print "before assignment",x_0,x
#    x  = x_0
#    y  = y_0
#    z  = z_0
#    vx = vx_0
#    vy = vy_0
#    vz = vz_0
#    print "after assignement",x_0,x
#
#print "during the header",x_0,x

# Opening Files
forces_total_file = open('output/Forces_total.csv','w')
forces_total_file.write("i value,Gravity X,Gravity Y,Gravity Z,Lift X,Lift Y,Lift Z,Drag X,Drag Y,Drag Z,Magnus X,Magnus Y,Magnus Z \n")
forces_total_file.close()

Moments_total_file = open('output/Moments_total.csv','w')
Moments_total_file.write("i value, alpha acceleration,phi acceleration, omega acceleration \n")
Moments_total_file.close()