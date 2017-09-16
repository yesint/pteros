.. pteros documentation master file, created by
   sphinx-quickstart on Thu Sep 25 11:57:14 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Python bindings for Pteros 
**************************

.. toctree::
	:maxdepth: 2
	
	system
	selection

General notes about the bindings
================================

Python bindings follow the C++ API of Pteros as close as possible, however they are desinged to look 'pythonic' rather than be a literal translation of C++ syntax. There are several major differences between C++ and Python APIs:

#. 
	The methods, which take output reference arguments in C++ return multiple values in Python instead. For example:
	
	.. code-block:: c++

		// C++
		Vector3f min,max;
		Selection.minmax(min,max);

	.. code-block:: python
	
		# Python
		(min,max) = Selection.minmax()

#. 
	All vector and matrix objects, which are returned by ``Pteros`` methods, are represented as *float* ``numpy`` arrays in Python instead of ``Eigen`` matrices. In addition in Python the vectors could be printed as rows without transposing them.

#. 
	In C++ the preferable way of accessing atom properties is using inlined accessor methods. In Python the preferable way is using indexing operator:

	.. code-block:: c++

		// C++
		Selection sel(sys,"name CA");
		cout << sel.Name(0) << endl;
		cout << sel.Resname(0) << endl;
		cout << sel.XYZ(0).transpose() << endl; // need to transpose to print vector as a row
	
	.. code-block:: python
		
		# Python
		sel = sys.select('name CA')
		print sel[0].name
		print sel[0].resname
		print sel[0].xyz # no need to transpose

 	The syntax ``sys.select(...)`` looks more natural in Python, but you can also use C++-like syntax ``Selection sys(...)`` if you like.

#. 
	There are C++-style accessor methods in Python, but they are prepended by *get* or *set* prefixes:
	
	.. code-block:: c++

		// C++		
		Selection sel(sys,"name CA");
		sel.Name(0) = "AH";          // write access
		cout << sel.Name(0) << endl; // read access
	
	.. code-block:: python
		
		# Python
		sel = sys.select("name CA")
		sel.setName(0,'AH')    # write access
		print sel.getName(0)   # read access

	Please note that there is no such thing as *l-value mathods* in Python, so you can't assing to accessor method as in C++.

Python bindings for the methods of Pteros classes could be either "fast" or "slow". "Fast" bindings map Pteros data structures directly to Python without any intermediate convertion thus their performance is almost the same as their C++ counterparts. "Slow" bindings use internal data convertion, which leads to significant performance penalty. Slow bindings are marked like this:

.. warning:: This method is slow due to...

You will probably never notice performance degradation in the simple scripts, but the usage of slow methods in the inner loops of complex algorithms may lead to unpleasant surprises. In such cases consider porting your algorithm to C++.


Using the bindings
==================

.. module:: pteros

Python bindings for Pteros are supplied as a single compiled extension library. Vectors and matrices in the bindings are represented as ``numpy`` objects, so most of the time you'll also need ``numpy`` module. In order to use the bindings just import ``pteros`` module and ``numpy``::

	from pteros import *
	import numpy as np


Classes
=======

* :class:`Atom`
* :class:`System`
* :class:`Selection`
* :class:`Periodic_box`
* :class:`Grid_searcher`