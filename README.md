# Landscape
Source codes for "Steering ecological-evolutionary dynamics to improve  artificial selection of microbial communities"
Directories in Codes are ordered according to figures.

## Figure 2:
* B: /Codes/Figure2/main_Figure2B_QuiverAttractor.m
* C: /Codes/Figure2/main_Figure2C_LandscapePlot.m
* E: See instructions for Figure 5C(ii) below


## Figure 3:
* A: Follow /Examples/Example1 to run the simulation codes in /Codes/Figure3/simu and generate the figure.
* B: The landscape and attractor are the same as in Figure 2(B, C).

## Figure 5:
* A(ii): /Codes/Figure2/main_Figure5Aii_Landscape2D.m
* B(ii, iv): /Codes/Figure5/B/ii_and_iv
 * Run main_repeat_simu.m to repeat the cycle 100 times.
 * To plot B(ii), run main_Figure5B_landscape.m.
 * To plot B(iv), run main_Figure5B_boxplot.m.
* B(iii, v): /Figure5/B/iii_and_v
 * Run main_repeat_simu.m to repeat the cycle 100 times.
 * To plot B(ii), run main_Figure5B_landscape.m.
 * To plot B(iv), run main_Figure5B_boxplot.m.
* C(ii): /Figure5/C/ii
 * Run main_heritability_simu_2cycle.m to simulate community functions of parents and offspring communities.
 * To plot C(ii), run main_Figure5C.m.
 * To plot Figure 2E, run main_Figure2E.m.
 * To plot Figure S4, run main_FigureS4.m.
* C(iii): /Figure5/C/iii
 * Run main_heritability_simu_2cycle_spike.m to simulate community functions of parents and offspring communities.
 * To plot C(iii), run main_Figure5C.m.

## Figure 6:
To generate the figures, follow /Examples/Example1 to run the simulation codes in /Codes/Figure6/simu with the following parameters:
  * A: C = 1000, comm_type_num = 2, Pn_sig = 0, spike_frac = 0;
  * B: C = 1000, comm_type_num = 2, Pn_sig = 0, spike_frac = 0.3;
  * C: C = 2000, comm_type_num = 10, Pn_sig = 0, spike_frac = 0;
  * D: C = 2000, comm_type_num = 10, Pn_sig = 0, spike_frac = 0.3;
  * E: C = 1000, comm_type_num = 2, Pn_sig = 100, spike_frac = 0;
  * F: C = 2000, comm_type_num = 2, Pn_sig = 100, spike_frac = 0.3;
  * G: C = 2000, comm_type_num = 10, Pn_sig = 100, spike_frac = 0;
  * H: C = 2000, comm_type_num = 10, Pn_sig = 100, spike_frac = 0.3;

## Figure 7:
To generate the figures, follow /Examples/Example2 to run the simulation codes in /Codes/Figure7/simu with the following parameters:
* A-B: check_period = C*10
* C-E: check_period = 100.

## Figure S2:
To generate the figures, follow /Examples/Example3 to run the simulation codes in /Codes/FigureS2_simu with the following parameters:
* top panel: end_cycle = 100, multiplier = 10;
* bottom panel: end_cycle = 100, multiplier = 100;

## Figure S3:
* A: /Codes/FigureS3/A/main_PCompare.m
* B: /Codes/FigureS3/B/main_PCompare.m

## Figure S4:
/Codes/Figure5/C/ii/
Run main_FigureS4.m after running main_heritability_simu_2cycle.m.

## Figure S5:
Bottom panel: Follow /Examples/Example1 to run the simulation codes in Codes/FigureS5/simu and generate the figure.

## Figure S6:
Follow instructions at the beginning of Codes/FigureS6/main_FigureS6.m to generate the figures.

## Figure S7:
Column 1 and 2: follow instructions at the beginning of Codes/FigureS7/main_Columns1and2.m to generate the figures.
Column 3: follow instructions at the beginning of Codes/FigureS7/main_Columns3.m to generate the figures.
Column 4: follow instructions at the beginning of Codes/FigureS7/main_Columns4.m to generate the figures.

## Figure S8: FigureS8/main_FigureS8.m

## Figure S9:
A: Follow /Examples/Example1 to run the simulation codes in Codes/FigureS9/simu and generate the figure.
B: Codes/Figure2/main_FigureS9B_Landscape2D.m

## Figure S10:
* A: To generate the figures, follow /Examples/Example1 to run the simulation codes in Codes/Figure6/simu with the following parameters:
  * T0 = 20, C = 1000, comm_type_num = 2, Pn_sig = 0, spike_frac = 0;
* B: Codes/FigureS10/BC/main_FigureS10B_landscape.m
* C: Codes/FigureS10/BC/
 Run main_heritability_simu_2cycle.m to simulate community functions of parents and offspring communities. Then plot C by running main_FigureS10C_heritbility.m.

## Figure S11:
Follow /Examples/Example1 to run the simulation codes in /FigureS11/simu and generate the figure. Before running main*.m, modify `spike_frac` accordingly. For example, to generate the figures of the second column, set spike_frac = 0.1. To reproduce the same dynamics, use the directory /10percentH.

## Figure S13:
To generate the figures, follow /Examples/Example2 to run the simulation codes in /Figure7/simu with the following parameters:
* A: HeriSwitch = int8(1), check_period = C*10. To reproduce the same dynamics, use the directory /Figure7/NoSpike.
* B: HeriSwitch = int8(1), check_period = 100. To reproduce the same dynamics, use the directory /Figure7/PeriodicUpdate.
* C: HeriSwitch = int8(-1), check_period = 100. To reproduce the same dynamics, use the directory /FigureS13/C.
* D: HeriSwitch = int8(0), check_period = 100. To reproduce the same dynamics, use the directory /FigureS13/D.

# Figure S14:
To generate the figures, follow /Examples/Example2 to run the simulation codes in /Figure7/simu with the following parameters:
* A: spike_initial = [**0.3** 0 0.6 -0.3 -0.6], check_period = C*10. To reproduce the same dynamics, use the directory /FigureS14/A.
* B: spike_initial = [**0.6** 0 0.3 -0.3 -0.6], check_period = C*10. To reproduce the same dynamics, use the directory /FigureS14/B.

## Figure S15:
To generate the figures, follow /Examples/Example2 to run the simulation codes in /Figure7/simu with the following parameters:
* A: HeriSwitch = int8(1), spike_clone_num = **1**. To reproduce the same dynamics, use the directory /FigureS15/A.
* B: HeriSwitch = int8(-1), spike_clone_num = **1**. To reproduce the same dynamics, use the directory /FigureS15/B.
* C: HeriSwitch = int8(0), spike_clone_num = **1**. To reproduce the same dynamics, use the directory /FigureS15/C.

## Figure S16:
To generate the figures, follow /Examples/Example2 to run the simulation codes in /Figure7/simu with the following parameters:
* A: HeriSwitch = int8(1), spike_clone_num = **2**. To reproduce the same dynamics, use the directory /FigureS16/A.
* B: HeriSwitch = int8(-1), spike_clone_num = **2**. To reproduce the same dynamics, use the directory /FigureS16/B.
* C: HeriSwitch = int8(0), spike_clone_num = **2**. To reproduce the same dynamics, use the directory /FigureS16/C.

## Figure S17:
To generate the figures, follow /Examples/Example2 to run the simulation codes in /Figure7/simu with the following parameters:
* A: HeriSwitch = int8(1), spike_clone_num = **10**. To reproduce the same dynamics, use the directory /FigureS16/A.
* B: HeriSwitch = int8(-1), spike_clone_num = **10**. To reproduce the same dynamics, use the directory /FigureS16/B.
* C: HeriSwitch = int8(0), spike_clone_num = **10**. To reproduce the same dynamics, use the directory /FigureS16/C.

## Figure S18:
To generate the figures, follow /Examples/Example2 to run the simulation codes in /FigureS18/simu. To reproduce the same dynamics, use the directory /FigureS18/Adaptive.

## Figure S19:
* A: /FigureS19/Landscape/A/main_LandscapePlot.m
* B: To generate the figures, follow /Examples/Example3 to run the simulation codes in /FigureS19/simu/B with end_cycle = 2e3 and multiplier = 100.
* C: To plot the landscape with the Newborn restrictor, run /FigureS19/Landscape/C/main_LandscapePlot.m.
To generate the evolutionary dynamics plots, follow /Examples/Example1 to run the simulation codes in /FigureS19/simu/CD with spike_frac = **0**. To reproduce the same dynamics, use the directory /FigureS19/C.
* D: To plot the landscape with the Newborn restrictor, run /FigureS19/Landscape/D/main_LandscapePlot.m.
To generate the evolutionary dynamics plots, follow /Examples/Example1 to run the simulation codes in /FigureS19/simu/CD with spike_frac = **0.3**. To reproduce the same dynamics, use the directory /FigureS19/D.

## Figure S20
* A: /FigureS20/A/main_LandscapePlot.m
* B: left: /FigureS20/BC_phiMt/main_phiM_t.m with spike_frac = 0 and phiM0 = [0.5 0.85].
    right: after running /FigureS20/D/main_heritability_simu.m, run /FigureS20/D/main_FigureS20B.m
* C: left: /FigureS20/BC_phiMt/main_phiM_t.m with spike_frac = 0.3 and phiM0 = [0.35 0.6].
  right: after running /FigureS20/E/main_heritability_simu.m, run /FigureS20/E/main_FigureS20C.m
* D: after running /FigureS20/D/main_heritability_simu.m, run /FigureS20/D/main_FigureS20D.m
* E: after running /FigureS20/E/main_heritability_simu.m, run /FigureS20/E/main_FigureS20E.m
* F: To generate the evolutionary dynamics plots, follow /Examples/Example1 to run the simulation codes in /Figure6/simu with the following parameters. To reproduce the same dynamics, use the directory /FigureS20/F.
 * C = 1000, T0 = **8.5**, comm_type_num = 10, Pn_sig = 0, spike_frac = **0**;
* G: To generate the evolutionary dynamics plots, follow /Examples/Example1 to run the simulation codes in /Figure6/simu with the following parameters. To reproduce the same dynamics, use the directory /FigureS20/G.
 * C = 1000, T0 = **8.5**, comm_type_num = 10, Pn_sig = 0, spike_frac = **0.3**;

## Figure S21
Follow README in /FigureS21.

## Figure S22
Follow README in /FigureS22.

## Figure S23
To generate the figures, follow /Examples/Example2 to run the simulation codes in /Figure7/simu with the following parameters:
* B: spike_initial = [0 0.3 -0.3], HeriSwitch = int8(1), spike_clone_num = 10, check_period = 100. To reproduce the same dynamics, use the directory /FigureS23/B.

## Figure S24
* No spiking, second cycle: run Figure5/C/ii/main_FigureS24top.m after running Figure5/C/ii/main_heritability_simu_2cycle.m.
* 30%-H spiking, second cycle: run Figure5/C/iii/main_FigureS24bottom.m after running Figure5/C/iii/main_heritability_simu_2cycle_spike.m.

## Figure S25
Figure S25 plots heritability check at Cycle 780 of the simulation plotted in black curves in FigureS17A.
* A, left panel: FigureS25/main_phiM_t.m with spike_frac = 0.6.
* B, rignt panel: FigureS25/main_phiM_t.m with spike_frac = 0.3.
To plot the figures in the 2nd and 3rd column, first run FigureS25/main_HM_SPNoise_HeriSwitch_Slope.m to reproduce the heritability check at Cycle 780. Then run FigureS25/main_Scatter_ParVSOff.m to generate the scatter plots.
