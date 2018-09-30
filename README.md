Multigrid Poisson Solvers
===============================================

Version 1 - December 2, 2017
by Matias di Martino <matiasdm@fing.edu.uy>


Introduction
------------

This is an implementation of several algorithms and criteria for Poisson
Editing.  The methods are detailed on the associated IPOL paper:

	"Multigrid Poisson Solvers"
	Matias di Martino and Gabriele Facciolo
	Image Processing On Line, 2016. DOI: XXX COMPLETE XXX
	http://dx.doi.org/ XXX COMPLETE XXX

Attention: The first time you execute MG.m, it will compute restriction and p
prolongation images and store them in a folder named "precompRP". This may take 
some time the first time you use the solver, but then, it allows the code to be 
faster as we don't need to calculate those matrix again. 


Files
-----

README.txt                  - This file.
LICENSE.txt                 - GNU AFFERO GENERAL PUBLIC LICENSE Version 3.
src/PoissonInpainting.m     - Solves poisson equation
src/lib/ConjugateGradint.m  - implementation of conjugate gradient
src/lib/createLaplacian...  - defines the laplacian operator
src/lib/GaussSeidel.m       - implementation of gaussSeidel method
src/lib/JacobiW.m           - weighted jacobi 
src/lib/MG.m                - Multigrid solver for the poisson eq.
examples.m                  - Examples of computation and use PoissonInpainting
images/*.png                - Images necessary for running the examples

The M-code files inside the folders "src/" and "src/lib" will be subjected to
the IPOL peer-review process.


Usage
-----

See the file "examples.m" for three examples of Poisson Image Editing using the
provided codes.


Portability
-----------

This implementation is intended to be compatible with ALL versions of Octave
and Matlab.

The only requirement is the Image Package for Octave or the Image Processing 
Toolkit for Matlab.  These requirements are only necessary for reading 

Any case of non-portability is considered a serious bug, and the authors would
like to be notified so that they can amend it.


Copyright and License
---------------------

Copyright (C) 2017,
 Matias Di Martino <matiasdm@fing.edu.uy>
 Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>

This is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

These files are distributed in the hope that they will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.



Thanks
------

The authors would be grateful to recieve any comment, especially about
portability issues, errors, bugs or strange results.
