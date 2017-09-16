System class
============

.. class:: System

	.. method:: __init__()
	
		Default constructor
		

	.. method:: __init__(sel_str)
	
		Constructor taking selection string as argument.
	
		:param str sel_str: Selection string
	

	.. method:: __init__(other_system)
	
		Constructor taking another system.
	
		:param System other_system: Other system to construct from.

		.. note:: Selections associated with the system are not copied!
		
	
	.. method:: num_atoms()
	
		:return: number of atoms in the system
		:rtype: int
		
	.. method:: num_frames()
	
		:return: number of frames in the system
		:rtype: int
	
	.. method:: select_all()
	
		:return: selection object for all atoms in the system
		:rtype: :class:`Selection`
		
	.. method:: select(str)
	
		:param string str: selection string
		:return: selection object for given selection string
		:rtype: :class:`Selection`
		
	.. method:: select(ind1,ind2)
	
		:param int ind1: first index
		:param int ind2: last index (inclusive)
		:return: selection object for given range of indexes
		:rtype: :class:`Selection`
		
	.. method:: select(l)
	
		:param sequence l: python list or tuple with selection indexes
		:return: selection object for given indexes
		:rtype: :class:`Selection`
		
	.. method:: load(fname[,b=0[,e=-1[,skip=-1[,on_frame=None]]]])
	
		Load any supported file format from disk (structure or trajectory)
		
		:param string fname: file to read
		:param int b: first frame to read
		:param int e: last frame to read
		:param int skip: read each n's frame
		:param on_frame: callable object to be called on each frame and passed current system as the first argument and the index of current frame as the second.
		:type on_frame: callable, f(System,int)
		
		Usage: ::
			
			# Define the callback
			def callback(sys,fr):
				print sys.num_atoms(),fr
				
			# Load structure file
			sys = System('struct.pdb')
			# Load trajectory with callback
			sys.load('traj.xtc',on_frame=callback)
			
		
	.. method:: frame_dup(fr)
	
		Duplicates given frame and adds it to the end of frame vector
	
	.. method:: frame_copy(fr1,fr2)
	
		Copy all frame data from fr1 to fr2
		
	.. method:: frame_delete([b=0[,e=-1]])
	
		Delete specified range of frames.
		If only *b* is supplied deletes all frames from *b* to the end.
		If only *e* is supplied deletes all frames from 0 to *e*.

		:param int b: first frame to delete
		:param int e: last frame to delete (inclusive)
    
	.. method:: getFrame_data(fr)    	
		
		:return: whole Frame structure from frame *fr*.
		:rtype: Frame object
		
	.. method:: setFrame_data(fr,data)
		
		:param Frame data: Frame object to set as frame *fr*
		
	.. method:: getBox(fr)    	
		
		:return: Periodic box from frame *fr*.
		:rtype: Periodic_box object
		
	.. method:: setBox(fr,data)
		
		:param Periodic_box data: box to set for frame *fr*
	
	.. method:: getTime(fr)    	
		
		:return: time stamp (in ps) from frame *fr*.
		:rtype: float
		
	.. method:: setTime(fr,t)
		
		:param float t: time stamp (in ps) to set for frame *fr*
		
	.. method:: getXYZ(ind,fr)
		
		:return: XYZ coordinates of atom *ind* from frame *fr*.
		:rtype: numpy.array(3)
		
	.. method:: setXYZ(ind,fr,coord)
		
		:param int ind: index of the atom
		:param int fr: frame to set
		:param coord: XYZ coordinates
		:type coord: array-like object of dimension 3 (list, tuple, numpy.array)
	
	.. method:: frame_append(frame)
		
		:param Frame frame: Frame object to add
		
	.. method:: assign_resindex()
		
		Assign unique resindexes. This is usually done automatically upon loading a structure from file.
	
	.. method:: atoms_dup(indexes)
		
		:param sequence indexes: list of tuple of atom indexes to duplicate
		
		.. warning:: This method is slow due to internal convertion between Python and C++ lists!
		
	.. method:: atoms_dup(atoms,coords)
		
		:param sequence atoms: list of tuple of atoms to add
		:param sequence coords: list of tuple of atom coordinates to add

		.. warning:: This method is slow due to internal convertion between Python and C++ lists!
		
	.. method:: wrap_all(fr[, dims=(1,1,1)])

		:param array_like(3) dims: dimensions to wrap
		
	.. method:: append(selection)
				append(system)

		:return: selection pointing to added atoms
		:rtype: Selection
		
	.. method:: dssp(filename)
	
		Determines secondary structure with DSSP algorithm and writes detailed report to file
		
	.. method:: dssp()
	
		Determines secondary structure with DSSP algorithm and returns it as a code string.
		The code is the same as in DSSP program:
		
		===========   ==============
		structure     code
		===========   ==============
		alphahelix	  "H"
		betabridge	  "B"
		strand		  "E"
		helix_3	      "G"
		helix_5	      "I"
		turn		  "T"
		bend		  "S"
		loop		  "\ " (space)
		===========   ==============
		
		:return: DSSP code string
                
	.. method:: sort_by_resindex()
	
		Sorts atoms by resindex arranging atoms with the same resindexes into contigous pieces. Could be called after atom additions or duplications.