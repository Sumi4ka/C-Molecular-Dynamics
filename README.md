# Name: C#: Molecular Dynamics
## Description
This code is written to solve problems of dynamics of complex molecular structures.
A feature of the solution is the consideration of the structure as an integral body 
with uniform translational and rotational characteristics.
This solution is implemented in C#.

This project consist of several classes:
- Atom and its child-classes (Example: Carbon)
- AtomRungeKutta4Coof that contain additional Runge Kutta coefficients for Atom
- Molecule that contain Atoms
- Like Atom Molecule has MoleculeRungeKutta4Coof
- Solver - main class for calculating
- Data - this class provide to work with Data:
  - In this catalog there is Data directory where there are two files with data

The program implements parallel programming.
And I use threadpool and hand-made synchronization.
