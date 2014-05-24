#!/usr/bin/python

from pteros import *

def cb(sys,fr):
    print 'callback reports:',fr,sys.num_atoms(),sys.num_frames()
    return True

s = System('2cfq.pdb')

s1 = System('2cfq.pdb')
print '>>> before append:',s1.num_atoms()
s1.append(System('2cfq.pdb'))
print '>>> after append:',s1.num_atoms()

print '>>> load with callback:'
s.load('2cfq.pdb',0,-1,0,cb)

print '>>> load without callback:'
s.load('2cfq.pdb')

b= s.getBox(0)
s.setBox(0,b)

coord = s.getXYZ(0,0)
print coord
coord[0] = 0
s.setXYZ(0,0,coord)
print s.getXYZ(0,0)

print s.dssp()

print 'before dup: ',s.num_atoms()
sel = Selection()
s.atoms_dup([1,2,3],sel)
print 'after dup: ',s.num_atoms()
print 'dup sel: ',sel.size()

sel = Selection(s,'name CA')
print '>>> sel.size:',sel.size()

