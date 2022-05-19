# pGFN-FF library

pGFN-FF Library
===============

This code contains a set of routines that build a library that is useful for generating
the parameters required to run a GFN-FF calculation or the periodic equivalent, pGFN-FF.

The starting point for the main part of the parameter generation is the original GFN-FF
force field of Spicher and Grimme, Angewandte Chemie Intl. Ed., 131, 11195 (2020) and
associated code within the XTB program (https://github.com/grimme-lab/xtb). 
This has been modified to allow for the changes proposed in the pGFN-FF method by 
Gale, LeBlanc, Spackman, Silvestri and Raiteri, J. Chem. Theory Comput., 17, 7827 (2021).

The code provide contains the routines to create the library and an example main program
that calls the library to generate the relevant parameters. 

License:

`pGFN-FF` is free software: you can redistribute it and/or modify it under the terms of the 
GNU Lesser General Public License as published by the Free Software Foundation, either 
version 3 of the License, or (at your option) any later version.

`pGFN-FF` is distributed in the hope that it will be useful, but without any warranty; 
without even the implied warranty of merchantability or fitness for a particular purpose. 
See the GNU Lesser General Public License for more details.
