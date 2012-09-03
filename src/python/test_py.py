import numpy
from pteros_py import *
fname = "/home/semen/work/Projects/pteros/svn/pteros/src/test/2lao.pdb"
s = System(fname)
print "Number of atoms: "+str(s.num_atoms())
s.load(fname,0)
print "Number of frames: "+str(s.num_frames())
p = s.Box(0)
print "Box:"
print p

print "Box Vectors & Angles:"
v = numpy.empty([3],dtype='f')
a = numpy.empty([3],dtype='f')
s.get_box_vectors_angles(0,v,a)
print v
print a

s.Time(1,0.02)
print s.Time(1)

s.frame_dup(0)
print s.num_frames()

print "--------------------------"
sel = Selection(s,"name CA")
print sel.num()
sel2 = sel
print sel2.num()
sel5 = Selection(sel)

sel.modify("name CA CB")
print sel.size()
print sel2.size()

sel3 = Selection(s,"name CA CB")

a = []
sel.each_residue(a)
print a[10].size()

print sel.get_system().num_atoms()
print sel.get_index()
print sel.get_chain()
print sel.get_unique_chain()
print sel.get_resid()
names = sel.get_name()
names[10] = "LOL"
print names
#sel.set_name(names)

print "Coordinates:"
arr = sel.get_xyz()
print arr
print "Testing set:"
arr[0,0] = 111
sel.set_xyz(arr)
print "Coordinates updated:"
print sel.get_xyz()

print "Average:"
print sel.get_average()

arr = sel.get_mass()
arr[0] = 112
sel.set_mass(arr)
print sel.get_mass()

trj = sel.get_traj(3)
print trj

print sel.center()

v1 = numpy.empty([3],dtype='f')
v2 = numpy.empty([3],dtype='f')

sel.rotate(1,0.2)
#sel.rotate(v1,v2)

print sel.rmsd(1)
print sel.rmsd(0,0)

print rmsd(sel,0,sel3,0)

sel.fit_trajectory()
print "After fit_trajectory:",sel.rmsd(1)

fit(sel3,sel)
print rmsd(sel,0,sel3,0)

print "Testing fit transform:"
sel.rotate(2,0.1)
tr = fit_transform(sel,sel3)
print tr
sel.apply_transform(tr)

sel.write("test.gro")

print "Testing inlined accessors:"
sel.set_frame(0)
print sel.X(10)
sel.setX(10,0,111)
print sel.X(10)

print sel.Resname(1)
sel.Resname(1,"TTT")
print sel.Resname(1)

a = sel.XYZ(4)
sel.setXYZ(4,a)
print sel.XYZ(4)

