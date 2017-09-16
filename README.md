# Pteros
Modern and fast molecular analysis and modeling library for C++ and Python

What is Pteros?
===============

Pteros is a C++ library for molecular modeling. It is designed to simplify the development of custom programs and scripts for molecular modeling, analysis of molecular dynamics trajectories and implementing new simulation and analysis algorithms. Pteros provides facilities, which are routinely used in all molecular analysis programs, namely input/output of popular file formats, powerful and flexible atom selections, geometry transformations, RMSD fitting and alignment, etc. Pteros also contains powerful facilities for parsing command-line arguments in custom programs and for running several analysis tasks in parallel, utilizing the power of modern multi-core processors.
Pteros supports writing analysis programs in either C++ or Python programming languages.

How to cite Pteros?
===================

Any work, which uses Pteros, should cite the following papers:

* Semen O. Yesylevskyy, "Pteros 2.0: Evolution of the fast parallel molecular analysis library for C++ and python", Journal of Computational Chemistry, 2015, 36(19), 1480–1488, doi: 10.1002/jcc.23943. (link)
* Semen O. Yesylevskyy, "Pteros: Fast and easy to use open-source C++ library for molecular analysis", Journal of Computational Chemistry, 2012, 33(19), 1632–1636. doi: 10.1002/jcc.22989. (link)

Is Pteros for you?
==================

Pteros library is for you if:

- You want to implement custom non-standard algorithms of molecular analysis.
- Your task is computationally expensive and potentially reusable.
- You want to run several "heavy" analysis tasks in parallel.
- You are not satisfied by the speed and memory consumption of the scripting languages embedded into popular molecular analysis programs, such as PyMol or VMD.
- You know C++ or don't mind learning this rather hard, but very powerful language.
- You know Python or want to learn it. Python scripts in Pteros are good for "throw-away" one-time scripts and serious reusable programs alike.

Pteros is not for you if:
=========================

- Your task requires extensive usage of molecular visualizer. Pteros doesn't have one currently.
- You have no programming skills at all, or you don't want to learn C++ or Python.

Features

- Reading/writing popular molecular structure and trajectory formats (PDB, GRO, MOL2, XTC, TRR, TPR, DCD, TNG).
- Very simple and expressive syntax for selecting groups of atoms similar to one used in VMD, but more powerfull.
- Selections can be manipulated in many different ways from simple translation to orientation by principal axes.
- Various properties of selections could be queried ranging from center of masses to radius of gyration.
- RMSD fitting and alignment.
- Full support for arbitrary periodic boxes - periodic selections and distances, wrapping/unwrapping, removing jumps over boundaries, computing closest images and shortest vectors, etc.
- Computing non-bonded interactions with any force field available in Gromacs format (GROMOS, CHARMM, AMBER, OPLS and more).
- Ability to work with very large trajectories, which does not fit into the memory.
- Asynchronous processing made easy. Many different analysis tasks could be run in parallel and simulataneously with trajectory reading.
- Very powerful and flexible syntax of the command-line options for custom analysis programs and scripts.
- Easy to use. API is very simple and intuitive.
- Easy to extend. Pteros is writen in high-level C++.
- Oriented to human beings (reserachers and students), not to programming gurus. The code is very well structured and documented.
- Bindings for the Python language.
- Plugin architecture, which makes writing C++ or Python code for asynchronous trajectory analysis trivial. You concentrate on the analysis algorithm itself, not on the technical details of managing files and trajectories.

