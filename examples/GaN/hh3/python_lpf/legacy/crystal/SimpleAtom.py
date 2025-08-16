

class Atom(object):

    def __init__(self, Z=None, symbol=None, mass=None):
        self.__dict__['Z'] = Z
        self.__dict__['symbol'] = symbol
        self.__dict__['mass'] = mass        
        return
