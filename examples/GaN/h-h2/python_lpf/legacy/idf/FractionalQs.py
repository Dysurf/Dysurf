#!/usr/bin/python

import numpy
from struct import pack,unpack,calcsize

version=1

intSize = calcsize('<i')
dubSize = calcsize('<d')
strSize = calcsize('<s')

def write(FQ,filename='FractionalQs',comment=''):
  """Takes numpy fractional Q with shape (N_q,D) and writes to binary file."""
  f=open(filename,'w')
  f.write(pack('<64s','FractionalQs'))
  f.write(pack('<i',version))
  f.write(pack('<1024s',comment))
  f.write(pack('<i',FQ.shape[1]))
  f.write(pack('<i',FQ.shape[0]))
  FQ = tuple( FQ.reshape(-1) )
  f.write(pack('<%id' % len(FQ),*FQ))
  return

def read(filename='FractionalQs'):
  """Takes filename, returns a tuple with information and fractional Q as a \n
     numpy."""
  f=open(filename,'r').read()
  i = 0
  filetype,= unpack('<64s',f[i:i+64*strSize])          ; i += 64*strSize
  version, = unpack('<i',f[i:i+intSize])               ; i += intSize
  comment, = unpack('<1024s',f[i:i+1024*strSize])      ; i += 1024*strSize
  D,N_q    = unpack('<2i',f[i:i+2*intSize])            ; i += 2*intSize
  FQ       = unpack('<%id' % (N_q*D),f[i:])
  FQ = numpy.array(FQ)
  FQ.shape = (N_q,D)
  return (filetype.strip('\x00'),version,comment.strip('\x00')),FQ

