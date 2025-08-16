#!/usr/bin/python

import numpy
from numpy import array
from struct import pack,unpack,calcsize

version=1

intSize = calcsize('<i')
dubSize = calcsize('<d')
strSize = calcsize('<s')

def write(Sq,filename='Sq',comment=''):
  """Takes numpy Sq in with shape (N_q) and writes to binary file."""
  f=open(filename,'w')
  f.write(pack('<64s','Sq'))
  f.write(pack('<i',version))
  f.write(pack('<1024s',comment))
  f.write(pack('<i',Sq.shape[0]))
  Sq = tuple( Sq.reshape(-1) )
  f.write(pack('<%id' % len(Sq),*Sq))
  return

def read(filename='Sq'):
  """Takes filename, returns a tuple with information and Sq as a numpy."""
  f=open(filename,'r').read()
  i = 0
  filetype, = unpack('<64s',f[i:i+64*strSize])          ; i += 64*strSize
  version,  = unpack('<i',f[i:i+intSize])               ; i += intSize
  comment,  = unpack('<1024s',f[i:i+1024*strSize])      ; i += 1024*strSize
  N_q,      = unpack('<i',f[i:i+1*intSize])             ; i += 1*intSize
  Sq       = unpack('<%id' % (N_q),f[i:])
  Sq = numpy.array(Sq)
  Sq.shape = (N_q)
  return (filetype.strip('\x00'),version,comment.strip('\x00')),Sq

