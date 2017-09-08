from _pteros11 import *


# Splitting selection by callback
def _split_by_callback(self,cb):
    m = {}
    parts = []
    for i in range(0,self.size()):
        ind = cb(self,i)
        if ind in m:
            m[ind].append( self[i].index )
        else:
            m[ind] = [ self[i].index ]
    for el in m:
        parts.append( Selection(self.get_system(),m[el]) )
    return parts

Selection.split = _split_by_callback
