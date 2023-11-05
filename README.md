# lightcurves

## use inputPhotons.py to generate afterglow phothons. 
You need to set the observed frequency: nuObs; total number of photon generated: totPh; number density of proton: nISM in the script


## compile the the photon-scattering.cpp and use the compiled executable to diffuse the photons in dense environments.
e.g. if the compiled executable names a.out, use

./a.out [parallel thread number] [input file name] [scale hight of disk in cgs] [number density in the disk midplanet in cgs] [viewing angle in deg] [half open angle of viewing cone in deg] [total number of photons needed within viewing cone]
