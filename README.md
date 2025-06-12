# dateset_of_nonequilibrium_cellular_uptake
This repository contains the code and data used in "Competition between activity-induced heterogeneity and phase interface in nonequilibrium cellular uptake".
# Langevin_dynamics
The Langevin dynamics simulation is performed using galamost, which is a open-source simulation tools employing GPUs (Zhu, Y.-L., et. al, J. Comput. Chem. 2013, 34, 2197-2211.). The algorithms of the membrane is based on the publications:Reynwar, B. J. et. al, Nature 2007, 447, 461-464. We implemented a custom membrane force field and the Berendsen barostat algorithm into galamost software.
# Equilibrium theory
The equilibrium membrane deformation was solved using the variational method implemented via MATLAB's bvp4c solver as implemented in bvp_2phase.m.
# Nonequilibrium theory
The evolution of the membrane under non-equilibrium conditions can be computed by running active_flow.m.
