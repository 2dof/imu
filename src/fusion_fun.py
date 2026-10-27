# Fusion Libraty for 6 dof IMU
# algorithms ported from 
#  1. https://github.com/xioTechnologies/Open-Source-AHRS-With-x-IMU.git    - only AHRS fusion update for 6DOF IMU
#  To DO:
#  2. Acceleration correction based on appendix A - Case I from 
#  " Summary of transformation equations And Equations of motions ued in free-flight
#     and wind-tunnel data reduction and analysis 
# Released under the MIT License (MIT)
# Copyright (c) 2025 Lukasz Szydlowski 

# Ver. 0.1 


from math import sqrt, sin, cos ,atan2, asin, degrees, radians, pi
 
#======================================================
def upadeIMU(q,a,g,BETA,Ts):
    # ported doem   MadgwickAHRS.c  
    # q =[q1,q2,q3,q4],              # quaternion vector 
    # a = [ax,ay,az],      # Acceleration [m/s^2
    # g [gx,gy,gz]         # Gyro  rad/s
     
    ax,ay,az=a
    gx,gy,gz=g
    q0,q1,q2,q3 =q 
    
    #// Rate of change of quaternion from gyroscope
    qDot1 = 0.5 * (-q1 * gx - q2 * gy - q3 * gz);
    qDot2 = 0.5 * (q0 * gx + q2 * gz - q3 * gy);
    qDot3 = 0.5 * (q0 * gy - q1 * gz + q3 * gx);
    qDot4 = 0.5 * (q0 * gz + q1 * gy - q2 * gx);
    
    #------ Accelerometer normalisaton
    anorm=sqrt(ax**2+ay**2+az**2)
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
    
    snorm=sqrt(s0**2+s1**2+s2**2+s3**2)

    #http://en.wikipedia.org/wiki/Fast_inverse_square_root
    inv_snorm=1/snorm

    s0 *= inv_snorm
    s1 *= inv_snorm
    s2 *= inv_snorm
    s3 *= inv_snorm
 
    #Apply feedback step
    qDot1 -= BETA * s0
    qDot2 -= BETA * s1
    qDot3 -= BETA * s2
    qDot4 -= BETA * s3
    #// Integrate rate of change of quaternion to yield quaternion
    q0 += qDot1 * Ts
    q1 += qDot2 * Ts
    q2 += qDot3 * Ts
    q3 += qDot4 * Ts
    
    qnorm=sqrt(q0**2+q1**2+q2**2+q3**2)
    inv_qnorm=1/qnorm
    
    q0 *= inv_qnorm
    q1 *= inv_qnorm
    q2 *= inv_qnorm
    q3 *= inv_qnorm

    return [q0,q1,q2,q3] 


def quaternion2euler(q):

    #return  pitch, roll, yaw  in radians
    # so calculate degress muliply by (180.0f / PI)
    # quaternion conjugate done in equations 
    a=180/pi
    pitch =-a* asin(2*q[1]*q[3]- 2*q[0]*q[2]); 
    roll  = a* atan2((2*q[2]*q[3] + 2*q[0]*q[1]), (2*q[0]*q[0] + 2*q[3]*q[3] -1))
    yaw   = a* atan2((2*q[1]*q[2] + 2*q[0]*q[3]), (2*q[0]*q[0]+ 2*q[1]*q[1] -1))
 
    return roll,pitch,yaw
def euler2quaternion(roll, pitch, yaw):
      """
      Convert an Euler angle to a quaternion.
      roll = roration of axis x
      pitch = roration of axis y
      yaw  = rotationn of axis z
      """
      r=roll/2
      p=pitch/2
      y=yaw/2
      qx = sin(r) * cos(p) * cos(y) - cos(r) * sin(p) * sin(y)
      qy = cos(r) * sin(p) * cos(y) + sin(r) * cos(p) * sin(y)
      qz = cos(r) * cos(p) * sin(y) - sin(r) * sin(p) * cos(y)
      qw = cos(r) * cos(p) * cos(y) + sin(r) * sin(p) * sin(y)
     
      return [qw,qx, qy, qz]
    
def quaternion_conj(q):
    return  [q[0], -q[1], -q[2],  -q[3]]


def calculate_tilt_angles(ax, ay, az): 
    """
    Calculate pitch and roll angles from acceleration.
    """
    roll = atan2(ay, sqrt(ax**2 + az**2)) * 180.0 / pi   # rotation of x
    pitch = atan2(-ax, az) * 180.0 / pi                       # rotation of y 
    return roll,pitch




if __name__ == '__main__':

# testing functions for:
# roll  = 45 deg
# pitch = 0.0
# yaw   = 0.0

    Ts=0.01
    
    # Statis test 
    roll =44.0  # deg
    a_z = cos(roll*pi/180)*9.8067
    a_y = sin(roll*pi/180)*9.8067
    a_x = 0.0

    a=(a_x, a_y, a_z)

    g=(0.0, 0.0, 0.0 ) #math.pi* (1 / 180.0) )

    BETA=0.05
     
    qi = euler2quaternion(roll*pi/180, 0.0, 0.0)

    #qi=[1.0, 0.0, 0.0, 0.0]
    

    qout = upadeIMU(qi,a,g,BETA,Ts)
    
 
     
    print('------------------------')
    print(quaternion2euler(quaternion_conj(qout)))
    #print(quaternion2euler(quaternion_conj(fus.q)))
    
    #------------------------------------------------
    #------------------------------------------------
    # imu bmi160 test with fusion
    from time import sleep_ms
    from machine import SoftI2C, Pin ,reset ,soft_reset
    import sys
    from machine import Pin, Timer
     
    import math
    from struct import unpack
    
    
    from bmi160 import BMI160
    
    #from ulab import numpy as np
    #print('ulab version', ulab.__version__)

    
    i2c = SoftI2C(scl=Pin(3), sda=Pin(2))      
    print('i2c adr:', i2c.scan())
    
    imu=BMI160(i2c)              # addr= 0x68  ( 108 dec)
    
    # test read
    print('raw ax, ay, az :',imu.read_raw_accel())
    print('raw wx, wy, wz :',imu.read_raw_gyro())
    # === read offset for acc 
    #data = imu.read_register(OFFSET_0, 3)
    acc_off = imu.readOffset_acc()

    # === read offset for gyto
    #data = imu.read_register(OFFSET_3, 3)
    gyro_off = imu.readOffset_gyro()
    
    
    imu.set_acc_off_en()
    print(' -------Auto calibration of ACC-------------')
    print('offset before: ',imu.readOffset_acc())
    print('acc before auto offset',imu.read_raw_accel())
    # 
    imu.foc_ax()
    sleep_ms(40)
    imu.foc_ay()
    sleep_ms(40)
    imu.foc_az()
    sleep_ms(40)
    
    imu.set_acc_off_en()
    if imu.read_bits(0x77, 6, 1): 
        print('auto acc offset: on')
    else: 
        print('auto acc offset : off')
    sleep_ms(100)
    print('offset: ',imu.readOffset_acc())
    print('acc calinrated',imu.read_raw_accel())
    
    # Gyro atuooffest
    imu.FOC_gyro()
    imu.set_gyro_off_en()
    
    # initialize  fusion 
    ws=3.14/180
    qi = euler2quaternion(0.0, 0.0, 0.0)
    
    ax,ay,az=imu.read_acceleration()
    wx,wy,wz=imu.read_gyro()
    
    #wx0,wy0,wz0= wx,wy,wz
    Ts=0.01
    print('-----------------')
    print('ax,ay,az,wx,wy,wz')
    #-----------
#     fus=Fusion()
#     fus.Ts=Ts
#     fus.q=qi    
    k=0

    for i in range(1000):
        k+=1
        
        #print(imu.read_gyro())
        #print(imu.read_acceleration())
        ax,ay,az = imu.read_acceleration()      # m/s^2
        wx,wy,wz=imu.read_gyro()                 # deg/sec
        
        wx *= 0.01744444 # convert to rad/sec
        wy *= 0.01744444
        wz *= 0.01744444 
        
        q = upadeIMU(qi,[ax,ay,az],[wx,wy,wz],BETA,Ts)
        qi= q
        
        #rx,ry,rz = quaternion2euler(quaternion_conj(q))
        rx,ry,rz = quaternion2euler(q)
        #print([rx,ry,rz])
        #fus.update_nomag([ax,ay,az],[wx,wy,wz])
        #rx1,ry1,rz1 = quaternion2euler(quaternion_conj(fus.q))
        
        
        #print('x:',rx,'y:',ry,'z:',rz)
        #print('z:',rz,'z1:',rz1)
        #print([ay,wy])
        if k>=10:
            k=0
            #print('ax',ax,'ay',ay)#,'az',az)
            rax,ray = calculate_tilt_angles(ax, ay, az)
            print('rx',rx,'rax',rax)
            #print('ry',ry,'ray',ray)
            #print(rx,ry)
#            print('x:',rx,'x1:',rx1,'wx',wx)
            #print('y:',ry,'y1:',ry1,'wy',wy)
            #print('z:',rz,'z1:',rz1,'wz',wz)
        sleep_ms(10)
        
    print(q)       
