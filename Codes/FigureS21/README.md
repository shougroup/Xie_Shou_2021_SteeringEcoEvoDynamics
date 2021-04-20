# Figure S21: simulating community selection of H-M community when the ecological interaction is mutualistic.
For the conditions in A and B, simulation is run 6 times for 6 independent replicates. To reproduce the results plotted in the figure, the random number generator state of each replicate is also provided. The simulation results are then processed and plotted.
## How to run
Simulation directories are organized as below:
* FigureS21/A
 *  FigureS21/A/1
     * FigureS21/A/1/C1
     * FigureS21/A/1/C2
     * ...
 *  FigureS21/A/2
     * FigureS21/A/2/C1
     * FigureS21/A/2/C2
     * ...
 * ...
* FigureS21/B
 *  FigureS21/B/1
     * FigureS21/B/1/C1
     * FigureS21/B/1/C2
     * ...
 *  FigureS21/B/2
     * FigureS21/B/2/C1
     * FigureS21/B/2/C2
     * ...
 * ...
where  /1 /2 ... are replicates, /C1 /C2 ... contains simulation results from Cycle 1, Cycle 2 ... To reproduce the same simulation results, /C1 of each simulation is provided according to the above structure.
1. To run replicate 1 of A, copy all the \*.m files in /FigureS21/simu into /FigureS21/A/1. Open the file \*ParaInti.m and modify the parameters accordingly. Set `s_check.min_check_int = 100` to check heritability every 100 cycles. Run the file main\*.m after modifying the number of cycles `C`.
2. After running all 6 replicates as described in step 1:
 * copy all \*.m files under /FigureS21/plot_simulations into /FigureS21/A.
 * Open main_PlotDynamics_Mutualistic.m, modify the number of cycles `C` and the range of directories `i` if necessary. If heritability check is performed, set `heri_flag = true`. Otherwise, set `heri_flag = false`.
 * Running main_PlotDynamics_Mutualistic.m will:
    * generate a directory /PlotData containing Data1.mat, Data2.mat, etc. Each \*.mat file contains evolutionary dynamics of each simulation replicate. If heritability checks are performed, /PlotData also contains directories Check1, Check 2 etc. Each /Check* directory contains results of heritability checks.
    * generate the figures in Figure 7C-E and Figure S13B.
4. Repeat 1-3 for /FigureS21/B where no spiking or heritability check is performed. In this case, set
 * `check_period = 6000` in \*ParaInti.m so that heritability is never checked.
 * `heri_flag = false` in main_PlotDynamics_Mutualistic.m.
5. To generate Figure S21C, run run main_strategy_compare.m.
## Description of output
 The codes save results from Cycle 1 into /C1, Cycle 2 into /C2, etc. Each /C\* folder contains:
 * comm_selected.mat: info of the selected communities.
 * distrng.mat: random number generator state before community reproduction.
 * H_isolates_in.mat: info of H clones of the spiking mix
 * M_isolates_in.mat: info of M clones of the spiking mix
 * /comm_all/P_all.mat: community function for all communities.
 * /comm_all/Pn.mat: measurement noise in community function, if noise is not zero.
 * /comm_all/newborns.mat: info of all Newborns in this cycle.  
