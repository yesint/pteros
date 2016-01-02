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
def _system_select(self,*args):
    if len(args)==0:
        # Select all
        sel = Selection(self,0,self.num_atoms()-1)
    else:
        sel = Selection(self)
        if len(args)==1:
            sel.modify(args[0])
        elif len(args)==2:
            sel.modify(args[0],args[1])
        else:
            raise RuntimeError("Wrong arguments for selection!")

    # keep system alive at the life time of selection
    sel._system = self    
    return sel

System.select = _system_select
System.__call__ = _system_select

#------------------------------------
# Wrapper for System.select_all()
#------------------------------------
def _system_select_all(self):
    sel = Selection(self,0,self.num_atoms()-1)
    # keep system alive at the life time of selection
    sel._system = self
    return sel

System.select_all = _system_select_all

#------------------------------------
# Wrapper for Selection.select() methods
#------------------------------------
def _selection_select(self,*args):
    sel = None
    if len(args)==1:
        sel = self._subselect(args[0])
    elif len(args)==2:
        sel = self._subselect(args[0],args[1])
    else:
        raise RuntimeError("Wrong arguments for subselection!")

    # keep system alive at the life time of selection
    sel._system = self.get_system()
    return sel

Selection.select = _selection_select
Selection.__call__ = _selection_select
