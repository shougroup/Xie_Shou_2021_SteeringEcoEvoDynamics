# Example1: simulating community selection of H-M community when only *f<sub>P</sub>*.
In most simulations, the same simulation is run 3 times for 3 independent replicates. To reproduce the results plotted in the figure, the random number generator state of each replicate is also provided. The simulation results are then processed and plotted.
## How to run
In this example, we use simulation codes from Figure6_simu to illustrate how to reproduce Figure 6A.
1. To run replicate 1, copy all the \*.m files in /Figure6/simu into /Figure6/A/1. Open the file main\*.m, modify the parameters according to the description. Press "Run" to run the simulation.
2. The codes save results from Cycle 1 into /C1, Cycle 2 into /C2, etc. Each /C\* folder contains:
  * comm_selected.mat: info of the selected communities.
  * distrng.mat: random number generator state before community reproduction.
  * /comm_all/P_all.mat: community function for all communities.
  * /comm_all/Pn.mat: measurement noise in community function, if noise is not zero.
  * /comm_all/newborns.mat: info of all Newborns in this cycle.  
3. After running all 3 replicates as described in step 1, copy all \*.m files under /Example1 into /Figure6/A. Open main_PlotDynamics.m, modify the number of cycles `C` and the range of `i` and hit "Run". This will generate a directory /PlotData containing Data1.mat, Data2.mat, etc. Each \*.mat file contains the info of each simulation replicate. It will also generate the two figures in Figure 6A.
