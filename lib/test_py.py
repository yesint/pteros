import numpy
from pteros_py import *

fname = "/home/semen/work/Projects/pteros/pteros/trunk/src/test/data/2lao.gro"
print ">>> Creating system from file"
s = System(fname)
print "Number of atoms: "+str(s.num_atoms())
s.load(fname)
print "Number of frames: "+str(s.num_frames())
print ">>> Get box"
p = s.getBox(0)
print "Box:"
print p
print ">>> Modify box[0,0] to 100 and read it back:"
p[0,0] = 100
print p
s.setBox(p,0)
print s.getBox(0)


print ">>> Box Vectors & Angles:"
v = numpy.empty([3],dtype='f')
a = numpy.empty([3],dtype='f')
s.get_box_vectors_angles(0,v,a)
print v
print a

print ">>> Set time of frame 1 to 0.02"
s.setTime(1,0.02)
print ">>> Get it back:"
print s.getTime(1)

print ">>> Duplicate frame 0 - number of frames now:"
s.frame_dup(0)
print s.num_frames()


print "--------------------------"
print ">>> Create selection 'name CA'"
sel = Selection(s,"name CA")
print "Size:", sel.num()
print ">>> Copy constructing another selection"
sel2 = sel
print "Size:", sel2.num()
print ">>> Assignment constructor"
sel5 = Selection(sel)
print "Size: ", sel5.size()

print ">>> Modify selection text to 'name CA CB'"
sel.modify("name CA CB")
print "Size now:", sel.size()
print "While copy is of size", sel5.size()
print "But sel2 is just a reference to sel1, thus size:",sel2.size()

print ">>> Make interval selection 10-12"
sel3 = Selection(s,10,12)
print "Size:",sel3.size()
print ">>> Modify it to another interval 20-30"
sel3.modify(20,30)
print "Size:",sel3.size()

print ">>> Testing each_residue in sel"
a = []
sel.each_residue(a)
print "Size of residue 10:", a[10].size()

print ">>> Testing various gets (10 first itemsa re shown)"
print "get_system().num_atoms():", sel.get_system().num_atoms()
print "get_index():", sel.get_index()[0:10]
print "get_chain():", sel.get_chain()[0:10]
print "get_unique_chain():",sel.get_unique_chain()[0:10]
print "get_resid():",sel.get_resid()[0:10]
print ">>> Get names, change 5th name to 'LOL', set back and get again"
names = sel.get_name()
names[5] = "LOL"
sel.set_name(names)
print sel.get_name()[0:10]

print ">>> Get coordinates:"
print sel.size()
arr = sel.get_xyz()
#for i in range(0,arr.size):
#    print i,arr[i]

print "Change (0,0) to 111 and show updated:"
arr[0,0] = 111
sel.set_xyz(arr)
print sel.get_xyz()

print ">>> Average:"
print sel.get_average()

print ">>> Getting mass"
arr = sel.get_mass()
print arr[0:10]
print ">>> Set first mass to 112 and get back"
arr[0] = 112
sel.set_mass(arr)
print sel.get_mass()[0:10]

print ">>> Getting trajectory for atom 3"
trj = sel.get_traj(3)
print trj

print ">>> Getting center of selection"
print sel.center()

print ">>> Rotate around x by 0.2"
#v1 = numpy.empty([3],dtype='f')
#v2 = numpy.empty([3],dtype='f')
sel.rotate(1,0.2)
#sel.rotate(v1,v2)

print ">>> RMSD of rotated frame 0 with not rotated frame 1"
print sel.rmsd(1)
print sel.rmsd(0,0)

#print rmsd(sel,0,sel3,0)

sel.fit_trajectory()
print ">>> After fit_trajectory:",sel.rmsd(1)

#fit(sel3,sel)
#print rmsd(sel,0,sel3,0)

print ">>> Testing transforms"
s = System()
s.load(fname)
s.frame_dup(0)
sel1 = Selection(s,"name CA")
sel2 = Selection(s,"name CA")
sel2.set_frame(1)
sel2.rotate(1,0.3)
print rmsd(sel1,0,sel2,1)
t = fit_transform(sel1,sel2)
sel1.apply_transform(t)
print rmsd(sel1,0,sel2,1)

print ">>> Testing write"
sel1.write("test.gro")

print ">>>Testing inlined accessors:"
sel1.set_frame(0)
print sel1.getX(10)
sel1.setX(10,0,111)
print sel1.getX(10)

print sel1.Resname(1)
sel1.Resname(1,"TTT")
print sel1.Resname(1)

a = sel1.getXYZ(4)
sel1.setXYZ(4,a)
print sel1.getXYZ(4)

print "--------------------------"
print "Testing Frame"
print "--------------------------"

l = s.getFrame_data(0).get_coord()
print l[0],l[1]

#print ">>> Testing Options_tree wrapper"
#o = Options_tree()
#o.from_command_line("--name aaa --type bbb --opt1 --sub1 1 val2 3 --end-opt1 --name ccc --name ddd")
#print o.to_command_line()
#print o.get_value_string(".type")
#print o.get_value_double(".haha",111.0)
#print o.get_values_int(".opt1.sub1")
#print o.get_option(".opt1").get_values_int("sub1")
#print o.get_options("name")

#print ">>> Testing Trajectory_processor"
#opt_str = "--trajectory_group /home/semen/work/Projects/kornelyuk/Sasha/dimer_md/1/dimer_noPBC_1.xtc --range frame_range 0 100 --end-range --structure_file /home/semen/work/Projects/kornelyuk/Sasha/dimer_md/1/dimer_pdb2gmx.gro --end-trajectory_group --async true --log_interval 100"
#opt = Options_tree()
#opt.from_command_line(opt_str)

#tr = Trajectory_processor(opt)
