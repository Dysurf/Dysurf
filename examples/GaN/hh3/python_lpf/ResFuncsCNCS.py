import numpy as N

class GaussianResolution :
    def __init__(self) :
        self._name='GaussianResolutionConvolution'

    def sigmae12(self,E):
    	""" Energy resolution function for CNCS for Ei=12meV
        	fitted to polynomial of order 4
        	-5.91725715232495e-07	5.63726943406554e-05	0.00100733850164346
        	-0.0461244825662230	0.476046292508122
    		"""
    	a0 =  0.476046292508122
    	a1 =  -0.0461244825662230
    	a2 =  0.00100733850164346
    	a3 =  5.63726943406554e-05
    	a4 =  -5.91725715232495e-07
    	sig = (a0 + a1*E +a2*E**2 + a3*E**3 +a4*E**4)/2.35 #convert to unit to c
    	return sig

    def sigmae20(self,E):
    	""" Energy resolution function for CNCS for Ei=20meV
        	fitted to polynomial of order 4
        	-1.64952077361938e-07	2.58258569481579e-05	0.000791281838645456
        	-0.0594951491635698	1.02188553396572
        """
    	a0 =  1.02188553396572
    	a1 =  -0.0594951491635698
    	a2 =  0.000791281838645456
    	a3 =  2.58258569481579e-05
    	a4 =  -1.64952077361938e-07
    	sig = (a0 + a1*E +a2*E**2 + a3*E**3 +a4*E**4)/2.35
    	return sig

    def sigmae20x5(self,E):
    	""" Energy resolution function for CNCS for Ei=20meV
        	fitted to polynomial of order 4
        	-1.64952077361938e-07	2.58258569481579e-05	0.000791281838645456
        	-0.0594951491635698	1.02188553396572
        """
    	a0 =  1.02188553396572
    	a1 =  -0.0594951491635698
    	a2 =  0.000791281838645456
    	a3 =  2.58258569481579e-05
    	a4 =  -1.64952077361938e-07
    	sig = (a0 + a1*E +a2*E**2 + a3*E**3 +a4*E**4)/2.35*5
    	return sig

def Gauss12(evalues,qvalues,ecenter,qcenter,sigmae12,sigmaq):
#    Ne=1.0/N.sqrt(2*N.pi*sigmae12(ecenter)**2)
    Ne=1.0/N.sqrt(2*N.pi*4*sigmae12(ecenter)**2)

    Nq=1.0/N.sqrt(2*N.pi*sigmaq**2)
    delta_q=qvalues-qcenter
    q_square=N.dot(delta_q,delta_q)
#    Res= Ne*N.exp(-(evalues-ecenter)**2/(2*sigmae12(ecenter)**2))*Nq*N.exp(-(q_square)/(2*sigmaq**2))
    Res= Ne*N.exp(-(evalues-ecenter)**2/(2*4*sigmae12(ecenter)**2))*Nq*N.exp(-(q_square)/(2*sigmaq**2))

    return Res
def Gauss20(evalues,qvalues,ecenter,qcenter,sigmae20,sigmaq):
    Ne=1.0/N.sqrt(2*N.pi*sigmae20(ecenter)**2)
    Nq=1.0/N.sqrt(2*N.pi*sigmaq**2)
    delta_q=qvalues-qcenter
    q_square=N.dot(delta_q,delta_q)
    Res= Ne*N.exp(-(evalues-ecenter)**2/(2*sigmae20(ecenter)**2))*Nq*N.exp(-(q_square)/(2*sigmaq**2))
    return Res

def BinFinder(binrange,value) :
    for i in range(len(binrange+1)) :
        if len(binrange) == 1:
            return 0
        if ((binrange[i] <= value) and (value <= binrange[i+1])) == True :
            #return (binrange[i+1]-binrange[i])/2.0
            return i
