# Liouville-von-Neumann-Matlab

[![View Liouville-von-Neumann-Matlab on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://uk.mathworks.com/matlabcentral/fileexchange/64703-liouville-von-neumann-matlab)

_Written by Marcin Konowalczyk and Gabriel Moise. [Timmel Group](http://timmel.chem.ox.ac.uk). University of Oxford; 2017_

This script performs a representative Liouville von Neumann simulation by propagation of the density matrix. The quantum mechanical system used for simulation consists of three spins: electrons (A,B) and a nucleus (C). Only one of the electrons is coupled to the nucleus (A-C) with a hyperfine coupling specified by `hfc`. The system is also subject to an external magnetic fields specified by `B0`. The calculation runs for the time points specified by `T`.
The code is intended to be used to learn about the basics of spin chemistry, not as a tool for simulation. It is heavily commented. To use it one should go though it line-by-line to understand what it does.
