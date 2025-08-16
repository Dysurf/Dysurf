#!/usr/bin/env python
# -*- coding: utf-8 -*-

from numpy import float64, cross, zeros, matrix, linalg
import numpy as np

#==============================================================#
# This code is used to convert the k/q-points in converntional #
# unit cell basis to primitive unit cell basis and vice versa  #
#==============================================================#

__author__ = 'Changpeng Lin'
__version__ = '0.1'
__email__ = 'lincp17@mails.tsinghua.edu.cn'
__date__ = 'Nov 6, 2018'

ifile = './uc_vector.in'
fi = open(ifile,'r')

con = zeros((3,3))
pri = zeros((3,3))
read_con = False
read_pri = False
read_kp = False
while True:
    line = fi.readline()
    if len(line) == 0:
        break
    if line.find('con') != -1:
        read_con = True
        n = 0
        continue
    if read_con == True:
        n += 1
        temp = float64(line.split())
        con[n-1] = temp
        if n == 3:
            read_con = False
    if line.find('pri') != -1:
        read_pri = True
        n = 0
        continue
    if read_pri == True:
        n += 1
        temp = float64(line.split())
        pri[n-1] = temp
        if n == 3:
            read_pri = False
    if line.find('c2p') != -1:
        mode = 'c2p'
        read_kp = True
        n = -1
        continue
    elif line.find('p2c') != -1:
        mode = 'p2c'
        read_kp = True
        n = -1
        continue
    if read_kp == True:
        n += 1
        if n == 0:
            nkq = int(line.split()[0])
            kkqq = zeros((nkq,3))
            continue
        temp = float64(line.split())
        kkqq[n-1] = temp

convol = con[0].dot(cross(con[1],con[2]))
privol = pri[0].dot(cross(pri[1],pri[2]))

bc = zeros((3,3))
bp = zeros((3,3))

bc[0] = cross(con[1],con[2])/convol
bc[1] = cross(con[2],con[0])/convol
bc[2] = cross(con[0],con[1])/convol

bp[0] = cross(pri[1],pri[2])/privol
bp[1] = cross(pri[2],pri[0])/privol
bp[2] = cross(pri[0],pri[1])/privol

if mode == 'c2p':       # Convert reduced coordinates from conventional to primitive
    bc = matrix(bc)*2*np.pi
    bp = matrix(bp)*2*np.pi
    convert = bc*linalg.inv(bp)   # Convertional matrix
    for  i in range(nkq):
        kkqq[i] = kkqq[i]*convert
        print(kkqq[i])
if mode == 'p2c':       # Convert reduced coordinates from primitive to conventional
    bc = matrix(bc)
    bp = matrix(bp)
    convert = bp*linalg.inv(bc)   # Convertional matrix
    for  i in range(nkq):
        kkqq[i] = kkqq[i]*convert
        print(kkqq[i])




