#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 14:21:16 2026

@author: hp
"""

from matplotlib.pylab import *
import matplotlib.pyplot as plt
plt.rcParams['lines.linewidth'] = 0.5
plt.rcParams['lines.markersize'] = 2
plt.rcParams["savefig.facecolor"]='white'
plt.rcParams["figure.facecolor"]='white'
plt.rcParams["axes.facecolor"]='white'
plt.ion()
plt.close('all')
import numpy as np
import math



#==================Fusion ========================
def quaternion_multiply(p,q): 
    # p : quaternion as: array 1x4 : [0,q1,q2,q3,q4]np
    # 
    # q : array 1x4 : [0,q1,q2,q3,q4]
    y=np.array([0.,0.,0.,0.])
    y[0] = (p[0] * q[0]) - (p[1] * q[1]) - (p[2] * q[2]) - (p[3] * q[3]);
    y[1] = (p[0] * q[1]) + (p[1] * q[0]) + (p[2] * q[3]) - (p[3] * q[2]);
    y[2] = (p[0] * q[2]) - (p[1] * q[3]) + (p[2] * q[0]) + (p[3] * q[1]);
    y[3] = (p[0] * q[3]) + (p[1] * q[2]) - (p[2] * q[1]) + (p[3] * q[0]);

    return y
 
def quaternion_modulus(q):
    return np.sqrt(np.sum(q*q))

def quaternion_conj(q):
    return np.array([q[0], -q[1], -q[2],  -q[3]])
   
def quaternion_norm(q):
    
    z=np.sqrt(np.sum(q*q))
    
    if z==0:
        return q
    return  q/z

def quaternion2euler(q):
    #return  pitch, roll, yaw  in radians
    # so calculate degress muliply by (180.0f / PI)
    a=180/np.pi
    pitch =-a* np.arcsin(2*q[1]*q[3]+ 2*q[0]*q[2]); 
    roll  = a* np.atan2((2*q[2]*q[3] - 2*q[0]*q[1]), (2*q[0]*q[0] + 2*q[3]*q[3] -1))
    yaw   = a* np.atan2((2*q[1]*q[2] - 2*q[0]*q[3]), (2*q[0]*q[0]+ 2*q[1]*q[1] -1))
        
    
    return pitch,roll,yaw    # rot of y,x,z

def euler2quaternion(roll, pitch, yaw):
      """
      Convert an Euler angle to a quaternion.
      roll = roration of axis x
      pitch = roration of axis y
      yaw  = oration of axis x
      """
      qx = np.sin(roll/2) * np.cos(pitch/2) * np.cos(yaw/2) - np.cos(roll/2) * np.sin(pitch/2) * np.sin(yaw/2)
      qy = np.cos(roll/2) * np.sin(pitch/2) * np.cos(yaw/2) + np.sin(roll/2) * np.cos(pitch/2) * np.sin(yaw/2)
      qz = np.cos(roll/2) * np.cos(pitch/2) * np.sin(yaw/2) - np.sin(roll/2) * np.sin(pitch/2) * np.cos(yaw/2)
      qw = np.cos(roll/2) * np.cos(pitch/2) * np.cos(yaw/2) + np.sin(roll/2) * np.sin(pitch/2) * np.sin(yaw/2)
     
      return [qw,qx, qy, qz]


def q2Rot(q):
    q0, q1, q2, q3 = q
    return np.array([
        [2*(q0**2 + q1**2) - 1, 2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2)],
        [2*(q1*q2 + q0*q3), 2*(q0**2 + q2**2) - 1, 2*(q2*q3 - q0*q1)],
        [2*(q1*q3 - q0*q2), 2*(q2*q3 + q0*q1), 2*(q0**2 + q3**2) - 1]
    ]) 

# #
GYRO_MEAN_ERROR = np.pi* (5.0 / 180.0) #5 deg/s gyroscope measurement error (in rad/s)  *from paper*

BETA =np.sqrt(3.0/4.0) * GYRO_MEAN_ERROR   # //*from paper*


DELTA_T=0.01  #// 100Hz sampling frequency
# ========================================
           
def quaternion_to_euler_zyx(q):
    """
    Convert quaternion to Z-Y-X Euler angles.

    Parameters
    ----------
    q : Quaternion
        Quaternion (in pyquaternion format)

    Returns
    -------
    np.array
        Z-Y-X Euler angles
    """

    # roll : rotation around x-axis
    sinr_cosp = 2 * (q[0] * q[1] + q[2] * q[3])
    cosr_cosp = 1 - 2 * (q[1] * q[1] + q[2] * q[2])
    roll = math.atan2(sinr_cosp, cosr_cosp)

    # pitch : rotation around y-axis
    sinp = 2 * (q[0] * q[2] - q[3] * q[1])
    if math.fabs(sinp) >= 1:
        pitch = math.copysign(math.pi / 2, sinp)
    else:
        pitch = math.asin(sinp)

    # yaw : rotation around z-axis
    siny_cosp = 2 * (q[0] * q[3] + q[1] * q[2])
    cosy_cosp = 1 - 2 * (q[2] * q[2] + q[3] * q[3])
    yaw = math.atan2(siny_cosp, cosy_cosp)

    # Euler angles
    return np.array([
        math.degrees(roll), 
        math.degrees(pitch), 
        math.degrees(yaw)
    ])

 

def upadeIMU(q,a,g,BETA,Ts):
    # ported doem   MadgwickAHRS.c  
    #See: http://www.x-io.co.uk/node/8#open_source_ahrs_and_imu_algorithms
    
    ax,ay,az=a
    gx,gy,gz=g
    q0,q1,q2,q3 =q 
    
    #	// Rate of change of quaternion from gyroscope
    qDot1 = 0.5 * (-q1 * gx - q2 * gy - q3 * gz);
    qDot2 = 0.5 * (q0 * gx + q2 * gz - q3 * gy);
    qDot3 = 0.5 * (q0 * gy - q1 * gz + q3 * gx);
    qDot4 = 0.5 * (q0 * gz + q1 * gy - q2 * gx);
    
    #------ Accelerometer normalisaton
    
    anorm=np.sqrt(ax**2+ay**2+az**2)
    inv_anorm=1/anorm
    ax *= inv_anorm
    ay *= inv_anorm
    az *= inv_anorm
    
    
    _2q0 = 2.0 * q0;
    _2q1 = 2.0 * q1;
    _2q2 = 2.0 * q2;
    _2q3 = 2.0 * q3;
    _4q0 = 4.0 * q0;
    _4q1 = 4.0 * q1;
    _4q2 = 4.0 * q2;
    _8q1 = 8.0 * q1;
    _8q2 = 8.0 * q2;
    q0q0 = q0 * q0;
    q1q1 = q1 * q1;
    q2q2 = q2 * q2;
    q3q3 = q3 * q3;
    
    #// Gradient decent algorithm corrective step
    s0 = _4q0 * q2q2 + _2q2 * ax + _4q0 * q1q1 - _2q1 * ay
    s1 = _4q1 * q3q3 - _2q3 * ax + 4.0 * q0q0 * q1 - _2q0 * ay - _4q1 + _8q1 * q1q1 + _8q1 * q2q2 + _4q1 * az;
    s2 =  4.0 * q0q0 * q2 + _2q0 * ax + _4q2 * q3q3 - _2q3 * ay - _4q2 + _8q2 * q1q1 + _8q2 * q2q2 + _4q2 * az;
    s3 =  4.0 * q1q1 * q3 - _2q1 * ax + 4.0 * q2q2 * q3 - _2q2 * ay;
    
    snorm=np.sqrt(s0**2+s1**2+s2**2+s3**2)
    
    #http://en.wikipedia.org/wiki/Fast_inverse_square_root
    
    inv_snorm= np.float32((1/snorm))
    
    s0 *= inv_snorm
    s1 *= inv_snorm
    s2 *= inv_snorm
    s3 *= inv_snorm
    
    # -- update 
    #// Apply feedback step
    qDot1 -= BETA * s0
    qDot2 -= BETA * s1
    qDot3 -= BETA * s2
    qDot4 -= BETA * s3
    
    
    #// Integrate rate of change of quaternion to yield quaternion
    q0 += qDot1 * Ts
    q1 += qDot2 * Ts
    q2 += qDot3 * Ts
    q3 += qDot4 * Ts
    
    qnorm=np.sqrt(q0**2+q1**2+q2**2+q3**2)
    inv_qnorm=np.float32(1/qnorm)
    
    q0 *= inv_qnorm
    q1 *= inv_qnorm
    q2 *= inv_qnorm
    q3 *= inv_qnorm

    return np.array([q0,q1,q2,q3])

data_path='/home/hp/Pobrane/madgwick_algorithm_matlab/ExampleData.mat'

import scipy.io
mat = scipy.io.loadmat(data_path)

gyr1=mat['Gyroscope']
acc1=mat['Accelerometer']
 
# plt.figure(1)
# plt.plot(gyr1)
# plt.figure(2)
# plt.plot(acc1)



GYRO_MEAN_ERROR = np.pi* (5.0 / 180.0) #5 deg/s gyroscope measurement error (in rad/s)  *from paper*
BETA =np.sqrt(3.0/4.0) * GYRO_MEAN_ERROR   # //*from paper*





 
Ts=0.01
BETA=0.01


roll =45

qi=euler2quaternion(roll*np.pi/180, 0.0, 0.0)

a_z=np.cos(roll*np.pi/180)*9.806
a_y=np.sin(roll*np.pi/180)*9.806
a_x=0.0

a=(a_x, a_y, a_z)
g=(0.0, 0.0, 0) #math.pi* (1 / 180.0) )

rxyz =quaternion_to_euler_zyx(qi)

 
qout = upadeIMU(qi,a,g,BETA,Ts)


# Testing 
# Qaternion to euler:
# 45 degree rotaton of X axis )roll

q_x45=euler2quaternion(45*np.pi/180, 0.0, 0.0)

ry,rx,rz =quaternion2euler(quaternion_conj(q_x45))

print('ry,rx,rz = pitch,roll,yaw [degree]:',ry,rx,rz)

R = q2Rot(q_x45)

# Rref =np.array([[ 1.00000107,  0. ,  0.  ],
#        [ 0.        ,  0.70710851, -0.70710634],
#        [ 0.        ,  0.70710634,  0.70710851]])
 


#=====================================================
 
from scipy.integrate import odeint

def process(x,t):
    #pendulum  - return in rsdians
    
    l=1
    y = x[0]
    dydt = x[1]
    dy2dt2 = - 9.806/l * np.sin(y)
    return [dydt,dy2dt2]  





# testing on modelled penulum 
roll = 45 
 
ns=150
delta_t = Ts
y = np.zeros((ns+1,2))  # process variable


#initial condition
w_x=0.0
roll = 90

yi = [roll*np.pi/180, w_x]
y[0] = yi
    
for i in range(0,ns):

    y0=odeint(process,yi,[0,delta_t])
    
    yi = y0[-1]
    y[i+1] = yi 


plt.figure(1)
plt.subplot(211)
plt.plot(y[:,0]*180/np.pi)
plt.ylabel('deg')
plt.grid(True)
plt.subplot(212)
plt.plot(y[:,1])
plt.ylabel('[r/s]')
plt.grid(True) 

a_z=np.cos(y[:,0])*9.806
a_y=np.sin(y[:,0])*9.806

gyrx=y[:,1]         


#--simulation with pendulum --------------
qout  = np.zeros((ns+1,4))  # process variable


qi=euler2quaternion(roll*np.pi/180, 0.0, 0.0)   # 4 

qout[0]=qi


rx0  = np.zeros(ns+1)  # process variable

rx0[0]=roll     # wyhhylenie pocztowe

for i in range(0,ns):
    
    ax = 0.0
    ay = a_y[i]
    az = a_z[i]
    
    gx = gyrx[i]
    gy = 0.0 
    gz = 0.0 
    
    qout[i+1] =  upadeIMU(qi,[ax,ay,ax],[gx,gy,gz],0.01 , delta_t)
    
    qi=qout[i+1]
    
    ryy,rx0[i+1],rzz =quaternion2euler(quaternion_conj(qi))
    
    print(y[i,0]*180/np.pi,rx0[i+1])
    
    
figure(2)
plt.plot(y[:,0]*180/np.pi,'r')
plt.plot(rx0)
plt.grid(True)

 
 

import struct

 
def fast_inverse_sqrt(number):
    
    
    # Convert the input number to a 32-bit floating point value
    x = float(number)
    
    # Interpret the input number as a 32-bit signed integer
    i = struct.unpack('i', struct.pack('f', x))[0]
    
    # Use the magic number to estimate the square root
    i = 0x5F3759DF - (i >> 1)
    
    # Interpret the estimate as a 32-bit floating point value
    y = struct.unpack('f', struct.pack('i', i))[0]
    
    # Use Newton's method to refine the estimate
    y = y * (1.5 - 0.5 * x * y * y)
    
    
    # Return the reciprocal square root
    return y


#import time  

# a= 0.0010 
# print(1/np.sqrt(a),fast_inverse_sqrt(a))

# d=1/np.sqrt(a)-fast_inverse_sqrt(a)
# print(d , ' ----  [%]', 100*d/(1/np.sqrt(a)))

# tic = time.perf_counter()
# 1/np.sqrt(a)
# toc = time.perf_counter()
# print('timed1:',toc - tic)

# tic = time.perf_counter()
# fast_inverse_sqrt(a)
# toc = time.perf_counter()
# print('timed2:',toc - tic)


# przesunięcie położenia imu wzglętem środks masy


 
