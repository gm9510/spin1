# spin1
Simulations of spin-1 magnetic systems

In this project c/c++ and Open MP was use in order to simulate the physics of low dimensional magnetic systems. Further studies and conclusions upon these results are presented in the book: Magnetic properties of spinor boson gases in optical lattices: a Gutzwiller variational approach.

Using the Gutzwiller ansatz, the imaginary time evolution and the runge-kutta method we found the ground state of the physical system of interest, see the book Magnetic properties of spinor boson gases in optical lattices: a Gutzwiller variational approach.

The code GA_sp_1.5.C look for the ground state for a given pair of parameters D and th.

The code GA_sp_2.8.C look for the ground state for a given value of th and several values of D

The code GA_sp_3.7.C look for the ground state for several values of D and th, and parallelize using OpenMP the process to optimize the time.
Huge amounts of data are generated from this process. The data is processed using the codes Glue_sp.C to order the data and Measure the observables.
