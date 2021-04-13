# Example2: simulating community selection of H-M community when 6 phenotypes can be modified by mutations.
In most simulations, the same simulation is run 3 times for 3 independent replicates. To reproduce the results plotted in the figure, the random number generator state of each replicate is also provided. The simulation results are then processed and plotted.
## How to run
In this example, we use simulation codes from Figure7/simu to illustrate how to reproduce Figure 7B.
1. To run replicate 1, copy all the \*.m files in /Figure7/simu into /Figure7/PeriodicUpdate/1. Open the file main\*.m, modify the parameters according to the description. Press "Run" to run the simulation.
2. The codes save results from Cycle 1 into /C1, Cycle 2 into /C2, etc. Each /C\* folder contains:
  * comm_selected.mat: info of the selected communities.
  * distrng.mat: random number generator state before community reproduction.
  * H_isolates_in.mat: info of H clones of the spiking mix
  * M_isolates_in.mat: info of M clones of the spiking mix
  * /comm_all/P_all.mat: community function for all communities.
  * /comm_all/Pn.mat: measurement noise in community function, if noise is not zero.
  * /comm_all/newborns.mat: info of all Newborns in this cycle.  
3. After running all 3 replicates as described in step 1:
 * copy all \*.m files under /Example2/plot_simulations into /Figure7/PeriodicUpdate.
 * Open main_PlotDynamics_SP.m, modify the number of cycles `C` and the range of directories `i`. If heritability check is performed, set `heri_flag = true`. Otherwise, set `heri_flag = false`.
 * Running main_PlotDynamics_SP.m will:
    * generate a directory /PlotData containing Data1.mat, Data2.mat, etc. Each \*.mat file contains evolutionary dynamics of each simulation replicate. If heritability checks are performed, /PlotData also contains directories Check1, Check 2 etc. Each /Check* directory contains results of heritability checks.
    * generate the figures in Figure 7C-E and Figure S13B.
4. Repeat 1-3 for /Figure7/NoSpike where no spiking or heritability check is performed. In this case, set
 * `check_period = C*10` in /Figure7/simu/main\*.m so that heritability is never checked.
 * `heri_flag = false` in main_PlotDynamics_SP.m.
5. To generate Figure 7F, copy all \*.m files under /Example2/compare_simulations into /Figure7. /Figure7 is the parent directory of both /NoSpike and /PeriodicUpdate. To compare the simulation results in these two folders, run main_strategy_compare.m.
