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


