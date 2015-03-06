** What is MinProj? **

MinProj is a MATLAB package that computes exact projection constants and minimal projections in coordinate spaces and matrix spaces, as well as approximate projection constants and minimal projections in polynomial spaces. It relies on the external software  CVX, and also on Chebfun for polynomial spaces. They are both included in the Basc folder you have downloaded.


** Which MinProj version is this? **

This is version 1.0.

** How do I start using Basc? **

Open MATLAB and go to the directory MinProj_v1.0 (the one containing this README.txt file). From there, activate the software by typing

MinProj_setup

You should then be able to call the generic MinProj commands (i.e., MinProjCoor, MinProjMatr, and MinProjPoly) from any directory. For instance, you can go into the subdirectory called documentation and run any portion of the reproducible.m file that comes with the article
"Computation of minimal projections and extensions"
by S. Foucart.

** Does MinProj come with any guarantee of accuracy? **

Absolutely not. MinProj is an academic exercise at its core and it contains undoubtedly many bugs at this stage. I encourage you to send me bug reports at

simon.foucart@centraliens.net

Please also use this address if you want to contact me about theoretical improvement or algorithm development for MinProj.
