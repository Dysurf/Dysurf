import sys
import numpy as N
import math
def rotation(vector,lattice):
    """Rotation matrix between Bin box (ql,qt1,qt2) and cartesion (qx,qy,qz)"""

    tolerance=1e-6
    rotation=N.zeros((3,3))
    x=[1,0,0]
    y=[0,1,0]
    z=[0,0,1]

    a1=lattice[0,:]
    a2=lattice[1,:]
    a3=lattice[2,:]
    vol=N.dot(a1,N.cross(a2,a3))

    b1=N.cross(a2,a3)/vol
    b2=N.cross(a3,a1)/vol
    b3=N.cross(a1,a2)/vol

#    b1L=math.sqrt(sum(b1**2))
#    b2L=math.sqrt(sum(b2**2))
#    b3L=math.sqrt(sum(b3**2))
    b1L=1.0
    b2L=1.0
    b3L=1.0

    if abs(vector[0]) <tolerance and abs(vector[1]) < tolerance:  # [001]
        rotation[2][0]=1.0
        rotation[0][1]=1.0
        rotation[1][2]=1.0
    elif abs(vector[0]) <tolerance and abs(vector[2]) < tolerance:  #[010]
        rotation[1][0]=1.0
        rotation[0][2]=1.0
        rotation[2][1]=1.0
    elif abs(vector[1]) <tolerance and abs(vector[2]) < tolerance:  #[100]
        rotation[0][0]=1.0
        rotation[1][1]=1.0
        rotation[2][2]=1.0
    elif abs(vector[0])<tolerance and abs(vector[1])>tolerance and abs(vector[2])>tolerance: # [0kl]
        theta= math.atan(b3L*vector[2]/(b2L*vector[1]))
        rotation[0][0]=1.0
        rotation[1][1]=math.cos(theta)
        rotation[1][2]=-1.0*math.sin(theta)
        rotation[2][1]=math.sin(theta)
        rotation[2][2]=math.cos(theta)
        rotation=N.dot(rotation,[z,y,x])     # change longitudinal to z

    elif abs(vector[0])>tolerance and abs(vector[1])<tolerance and abs(vector[2])>tolerance: # [h0l]
        theta= math.atan(b3L*vector[2]/(b1L*vector[0]))
        rotation[1][1]=1.0
        rotation[0][0]=math.cos(theta)
        rotation[0][2]=math.sin(theta)
        rotation[2][0]=-1.0*math.sin(theta)
        rotation[2][2]=math.cos(theta)

    elif abs(vector[0])>tolerance and abs(vector[1])>tolerance and abs(vector[2])<tolerance: # [hk0]
        theta= math.atan(b2L*vector[1]/(b1L*vector[0]))
        rotation[2][2]=1.0
        rotation[0][0]=math.cos(theta)
        rotation[0][1]=-1.0*math.sin(theta)
        rotation[1][0]=math.sin(theta)
        rotation[1][1]=math.cos(theta)
#    elif abs(vector[0] -vector[1])< tolerance and abs(vector[1] - vector[2]) < tolerance and abs(vector[2] - vector[0]) <tolerance: # [111]
    else:
        sys.exit("Input hkl is wrong,  there should be at lest one zero in (hkl)")

    return rotation
