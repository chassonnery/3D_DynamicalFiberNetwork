The Fortran files can be compilated into a executable file **analysis.exe** using the command line :  
> `gfortan definitions.f90 fonctions.f90 lecture_ecriture.f90 jacobi.f90 analysis.f90 -o analysis.exe`

or 
> `gfortan definitions.f90 fonctions.f90 lecture_ecriture.f90 jacobi.f90 analysis_10runs.f90 -o analysis.exe`

The first alternative will ask the user to enter the path to a folder containing the result of *one simulation*. It will then compute the alignment rate, number of links and number of neighbours of each fiber at each time-step, as well as the average of these quantifiers over all fibers at each time-step.

The second alternative will ask the user to enter the path to a folder containing the result of *10 simulations*, each contained in a subfolder numbered from run0 to run9. It will then perform the analysis described above for each simulation, as well as compute the average quantifiers over 10 runs at each time-step. This only give meaningfull results if the 10 simulations have the same parameters.

Finally, the script xxx.py can be used to reproduce the graphics presented in the article `Fiber crosslinking drives the emergence of order in a 3D dynamical network model`.
