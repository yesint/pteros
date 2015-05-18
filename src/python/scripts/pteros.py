# This is main entry point for pteros module
from _pteros import *

#=========================================================================
# Bindings for System.select() are implemented on Python side because
# otherwise expressions like
# System('file.pdb').select('name CA').write('res.gro')
# do not work. The reason is that C++ destructor of System is called
# *before* Python calls write(). In C++ this never happens.
# To prevent this we inject a variable _system to returned selection,
# which holds a reference to the parent System until selection is alive.
# This is ugly, but solves the problem.
#=========================================================================

#------------------------------------
# Wrapper for System.select() methods
#------------------------------------
def _select(self,*args):
    sel = Selection(self)
    if len(args)==1:
        sel.modify(args[0]) # For string or index array
        # keep system alive at the life time of selection
        sel._system = self
        return sel
    elif len(args)==2:
        sel.modify(args[0],args[1]) # For pair of indexes
        # keep system alive at the life time of selection
        sel._system = self
        return sel

System.select = _select
System.__call__ = _select

#------------------------------------
# Wrapper for System.select_all()
#------------------------------------
def _select_all(self):
    sel = Selection(self,0,self.num_atoms()-1)
    # keep system alive at the life time of selection
    sel._system = self
    return sel

System.select_all = _select_all
