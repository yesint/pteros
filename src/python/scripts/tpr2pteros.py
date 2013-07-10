#!/usr/bin/python

import re
import sys
import re
import os
from subprocess import call

print "*** This is tpr2pteros.py script ***"

class Atom:
	pass

# Check if we have pttop file already
if os.path.isfile(sys.argv[1]+".pttop"):
	# Check modification time of pttop and tpr and decide if need to proceed
	tpr_time = os.path.getmtime(sys.argv[1])
	pttop_time = os.path.getmtime(sys.argv[1]+".pttop")
	if(tpr_time < pttop_time):
		print "PTTOP file '%s' is up to date. Exiting." % (sys.argv[1]+".pttop")
		sys.exit()

# Run gmxdump
print "Calling gmxdump on tpr file '%s'..." % sys.argv[1]
call("gmxdump -sys -s %s > _tpr_dump" % sys.argv[1], shell=True)

# Read the dump
f = open("_tpr_dump","r")

atoms = []
cg = []
excl = []
LJ = []
LJ14 = []
LJ14_pairs = []
bonded_pairs = []
fudgeQQ = 0.0
box = []

while True:
	line = f.readline()
	if not line: break
	if re.search(" atoms:",line):
		print "Reading atoms..."
		line = f.readline()
		m = re.search("(\d+)",line)
		natoms = int(m.group(1))
		print "natoms = %s" % natoms
		for i in range(0,natoms):
			line = f.readline()
			m = re.search("type=\s*(\d+),.+m=\s*(\S+),.+q=\s*(\S+),.+resind=\s*(\d+),",line)
			at = Atom()
			at.type = int(m.group(1))
			at.mass = float(m.group(2))
			at.charge = float(m.group(3))
			at.resind = int(m.group(4))
			atoms.append(at)
			
		# Now read atom names
		line = f.readline()
		for i in range(0,natoms):
			line = f.readline()
			m = re.search("\"(\S+)\"",line)
			atoms[i].name = m.group(1)
			
		# Now read atom typenames
		line = f.readline()
		for i in range(0,natoms):
			line = f.readline()
			m = re.search("name=\"(\S+)\",",line)
			atoms[i].typename = m.group(1)
			
		# Now read number of residues
		line = f.readline()
		m = re.search("(\d+)",line)
		nres = int(m.group(1))
		print "Number of residues: %s" % nres
		# Read residue names
		last = 0
		for i in range(0,nres):
			if i % 1000 ==0:
				print i
			line = f.readline()			
			m = re.search("name=\"(\S+)\",.+nr=(\d+)",line)
			resname = m.group(1)
			resid = int(m.group(2))
			# assign to atoms
			for ind in range(last,len(atoms)):
				#print i,atoms[ind].resind
				if atoms[ind].resind==i:
					atoms[ind].resid = resid
					atoms[ind].resname = resname					
				else:
					#print "--------"
					last = ind
					break
					
			#for at in atoms:
			#	if at.resind==i:
			#		at.resid = int(resid)
			#		at.resname = resname

		continue
			
	if re.search("cgs:",line):
		print "Reading charge groups..."
		# Charge groups
		line = f.readline()
		m = re.search("(\d+)",line)
		ncgs = int(m.group(1))
		print "There are %s charge groups" % ncgs
		for i in range(0,ncgs):
			line = f.readline()			
			m = re.search("(\d+)\.\.(\d+)",line)
			cg.append( (int(m.group(1)),int(m.group(2))) )

	if re.search("excls:",line):
		print "Reading exclusions..."
		# Exclusions
		line = f.readline()
		m = re.search("(\d+)",line)
		nexcl = int(m.group(1))
		line = f.readline()
		print "There are %s atoms with exclusions" % nexcl
		for i in range(0,nexcl):
			line = f.readline()
			m = re.search("\{(.+)",line) #Select everything starting from {
			# Split line
			l = re.split(",\s*|\}",m.groups(1)[0])
			# Check for continuations
			if not re.search("\}",line):
				line = f.readline()
				m = re.search("(.+)\}",line) #Select everything before }
				l += re.split(",\s*",m.groups(1)[0])
			# now transform to numbers
			excl.append( [int(x) for x in l if x!=''] )
			
	if re.search("idef",line):
		print "Reading interactions..."
		# Interaction definitions
		line = f.readline()
		m = re.search("(\d+)",line)
		atnr = int(m.group(1)) # Number of distinct non-bond atom types
		line = f.readline()
		m = re.search("(\d+)",line)
		ntypes = int(m.group(1)) # Number of interaction types
		# Next antr*atnr lines contain non-bond matrix of LJ_SR interactions
		# They are used for all pairs, which are not listed explicitly in exclusions		
		r = []
		for i in range(0,atnr*atnr):
			line = f.readline()
			m = re.search("c6=\s*(.+),\s+c12=\s*(.+)",line)
			c6 = float(m.group(1))
			c12 = float(m.group(2))
			# Here is how access is given in GROMACS between atoms with atomtypes A and B:
			# float c6=mtop->ffparams.iparams[B*atnr+A].lj.c6;
			# float c12=mtop->ffparams.iparams[B*atnr+A].lj.c12;
			r.append( (c6,c12) )
				
		for A in range(0,atnr):
			el = []
			for B in range(0,atnr):
				el.append(0)
			LJ.append(el)
				
		for A in range(0,atnr):
			for B in range(0,atnr):
				LJ[A][B] = r[B*atnr+A]
				LJ[B][A] = r[B*atnr+A]
		
		# Now read all other interactions
		# For now we only extract LJ14 from here
		for i in range(0,ntypes-atnr*atnr):
			line = f.readline()
			m = re.search("functype\[(.+)\]=(\S+),\s*\S+=\s*(\S+),\s*\S+=\s*(\S+),",line)			
			if not m:
				continue
			n = int(m.group(1))
			if m.group(2) == "LJ14":
				c6 = float(m.group(3))
				c12 = float(m.group(4))
				LJ14.append( (n,c6,c12) )
				
		# Read fudgeQQ:
		line = f.readline()
		m = re.search("=\s*(.+)",line)
		fudgeQQ = float(m.group(1))
		
	if re.search("LJ-14:",line):
		# Read LJ14 pairs
		line = f.readline()
		line = f.readline()
		while True:
			line = f.readline()
			m = re.search("type=(\d+)\s*\((\S+)\)\s+(\d+)\s+(\d+)",line)
			if not m:
				break
			if m.group(2) == "LJ14":
				LJ14_pairs.append( (int(m.group(3)), int(m.group(4)), int(m.group(1))) )
	
	# Read all bonded information
	if re.search("Bond:",line):
		# Read constraines
		line = f.readline()
		line = f.readline()
		while True:
			line = f.readline()
			m = re.search("type=(\d+)\s*\((\S+)\)\s+(\d+)\s+(\d+)",line)
			if not m:
				break
			if m.group(2) == "BONDS":
				bonded_pairs.append( (int(m.group(3)), int(m.group(4))) )
	
	if re.search("Constraint:",line):
		# Read constraines
		line = f.readline()
		line = f.readline()
		while True:
			line = f.readline()
			m = re.search("type=(\d+)\s*\((\S+)\)\s+(\d+)\s+(\d+)",line)
			if not m:
				break
			if m.group(2) == "CONSTR":
				bonded_pairs.append( (int(m.group(3)), int(m.group(4))) )
	
	if re.search("Settle:",line):
		# Read settles for water and convert them to bonds	
		line = f.readline()
		line = f.readline()
		while True:
			line = f.readline()
			m = re.search("type=(\d+)\s*\((\S+)\)\s+(\d+)\s+(\d+)\s+(\d+)",line)
			if not m:
				break
			if m.group(2) == "SETTLE":
				bonded_pairs.append( (int(m.group(3)), int(m.group(4))) )
				bonded_pairs.append( (int(m.group(3)), int(m.group(5))) )	
	
	if re.search("box \(",line):
		# Periodic box
		for i in range(0,3):
			line = f.readline()
			m = re.search("\{\s*(\S+),\s*(\S+),\s*(\S+)\}",line)
			box.append( float(m.group(1)) )
			box.append( float(m.group(2)) )
			box.append( float(m.group(3)) )
			

	if re.search("^x \(\d+x\d+\)",line):
		print "Reading coordinates..."
		# Coordinates
		for i in range(0,len(atoms)):
			line = f.readline()
			m = re.search("\{\s*(\S+),\s*(\S+),\s*(\S+)\}",line)
			x = float(m.group(1))
			y = float(m.group(2))
			z = float(m.group(3))
			atoms[i].coord = (x,y,z)

#-----------------------------------------
# Output
#-----------------------------------------
f = open(sys.argv[1]+".pttop","w")

f.write("# number of atoms\n")
f.write("%s\n" % len(atoms))
f.write("# N atomname typename atomtype resid resind resname mass charge x y z\n")
i = 0
for	at in atoms:
	f.write("%i %s %s %i %i %i %s %e %e %e %e %e\n" % (i,at.name, at.typename, at.type, at.resid, at.resind, at.resname, at.mass, at.charge, at.coord[0], at.coord[1], at.coord[2]) )
	i += 1

f.write("# number of bonds\n")
f.write("%i\n" % len(bonded_pairs))
f.write("# bonds (atom1,atom2)\n")
for b in bonded_pairs:
	f.write("%i %i\n" % (b[0],b[1]))
	
f.write("# Periodic box\n")	
n = 0
for i in range(0,3):
	for j in range(0,3):
		f.write("%e " % box[n])		
		n += 1
	f.write("\n")

f.write("# number of charge groups\n")
f.write("%s\n" % len(cg))
f.write("# charge groups (begin,end)\n")
for	c in cg:
	f.write("%i %i\n" % (c[0],c[1]) )


f.write("# number of atoms with exclusions\n")
f.write("%s\n" % len(excl))
f.write("# exclusions: atom --> excluded atoms\n")
i = 0
for	e in excl:
	f.write("%i\t" % i)
	for a in e:
		f.write("%i " % a)
	f.write("\n")
	i += 1

f.write("# size of LJ matrix\n")
f.write("%s\n" % len(LJ))
f.write("# LJ matrix of C6\n")
for	row in LJ:
	for el in row:
		f.write("%e\t" % el[0])
	f.write("\n")
f.write("# LJ matrix of C12\n")
for	row in LJ:
	for el in row:
		f.write("%e\t" % el[1])
	f.write("\n")

f.write("# fudgeQQ\n")
f.write("%s\n" % fudgeQQ)

f.write("# number of LJ distinct interaction types\n")
f.write("%s\n" % len(LJ14))
f.write("# LJ14 types: number, C6, C12\n")
table = {}
i = 0
for	e in LJ14:
	table[e[0]] = i
	f.write("%i\t%e\t%e\n" % (i,e[1],e[2]))
	i += 1

f.write("# number LJ14 atom pairs\n")
f.write("%s\n" % len(LJ14_pairs))
f.write("# LJ14 atom pairs: i, j, type_number\n")
for	e in LJ14_pairs:
	f.write("%i\t%i\t%i\n" % (e[0],e[1], table[e[2]]))

f.close()

# Remove temporary file if not asked to keep it
if not sys.argv[2]=="keep":
	os.remove("_tpr_dump")
