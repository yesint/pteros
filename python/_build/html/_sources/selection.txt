Selection class
===============

.. function:: rmsd(selection1, frame1, selection2, frame2)
	
	Computes RMSD of two selections of the same size for given frames
	
	:return: RMSD value
	:rtype: float
	
.. function:: fit(selection1, selection2)
	
	Fits *selection1* to *selection2*

.. function:: fit_transform(selection1, selection2)
	
	Computes fit transform for fitting *selection1* to *selection2*
	
	:return: fit transform matrix
	:rtype: 4x4 matrix

.. function:: non_bond_energy(selection1, selection2 [, cutoff=0.25 [, fr=-1 [, periodic=True]]])
	
	Computes non-bond interaction energy between two selections within given cut-off distance.
	If *fr* is not specified the current frame of *selection1* is used.

.. class:: Selection

	.. method:: __init__(system, sel_str)
	
		Create selection from selection string
		
	.. method:: __init__(system)
	
		Create empty selection bound to system
	
	.. method:: __init__(selection)
	
		Create selection from other selection
	
	.. method:: __init__(system, ind1, ind2)
	
		Create selection from the range of absolute indexes
	
	.. method:: __init__(system, seq)
	
		Create selection from the sequence of indexes
		
		:param sequence seq: list or tuple of indexes
		
	.. method:: size()
	
		Returns selection size

        .. method:: append(selection)

                Appends another selection to this one.

        .. method:: append(ind)

                Appends absolute index to selection.

                :param int ind: index to append

        .. method:: remove(sel)

                Removes all atoms of sel from current selection

                :param Selection sel: selection to remove

        .. method:: remove(ind)

                Removes absolute index from selection.

                :param int ind: index to remove

        .. method:: set_system(system)

                Sets new system for selection.

                .. warning:: This clears selection index and leaves it empty!
