# Landscape
Source codes for "Steering ecological-evolutionary dynamics to improve  artificial selection of microbial communities"
Directories in Codes are ordered according to figures.

## Figure 2:
* B: /Figure2/main_Figure2B_QuiverAttractor.m
* C: /Figure2/main_Figure2C_LandscapePlot.m
* D: See instructions for Figure 5C(ii) below


## Figure 3:
* A: Follow /Examples/Example1 to run the simulation codes in /Figure3 and generate the figure.
* B: The landscape and attractor are the same as in Figure 2(B, C).

## Figure 5:
* A(ii): /Figure2BC/main_Figure5Aii_Landscape2D.m
* B(ii, iv): /Figure5/B/ii_and_iv
 * Run main_repeat_simu.m to repeat the cycle 100 times.
 * To plot B(ii), run main_Figure5B_landscape.m.
 * To plot B(iv), run main_Figure5B_boxplot.m.
* B(iii, v): /Figure5/B/iii_and_v
 * Run main_repeat_simu.m to repeat the cycle 100 times.
 * To plot B(ii), run main_Figure5B_landscape.m.
 * To plot B(iv), run main_Figure5B_boxplot.m.
* C(ii): /Figure5/C/ii
 * Run main_heritability_simu.m to simulate community functions of parents and offspring communities.
 * To plot C(ii), run main_Figure5C.m.
 * To plot Figure 2D, run main_Figure2D.m.
 * To plot Figure S4, run main_FigureS4.m.
* C(ii): /Figure5/C/iii
 * Run main_heritability_simu.m to simulate community functions of parents and offspring communities.
 * To plot C(iii), run main_Figure5C.m.

## Figure 6:
To generate the figures, follow /Examples/Example1 to run the simulation codes in /Figure6_simu with the following parameters:
  * A: C = 1000, comm_type = 2, Pn_sig = 0, spike_frac = 0;
  * B: C = 1000, comm_type = 2, Pn_sig = 0, spike_frac = 0.3;
  * C: C = 2000, comm_type = 10, Pn_sig = 0, spike_frac = 0;
  * D: C = 2000, comm_type = 10, Pn_sig = 0, spike_frac = 0.3;
  * E: C = 1000, comm_type = 2, Pn_sig = 100, spike_frac = 0;
  * F: C = 2000, comm_type = 2, Pn_sig = 100, spike_frac = 0.3;
  * G: C = 2000, comm_type = 10, Pn_sig = 100, spike_frac = 0;
  * H: C = 2000, comm_type = 10, Pn_sig = 100, spike_frac = 0.3;

## Figure 7:
To generate the figures, follow /Examples/Example2 to run the simulation codes in /Figure6_simu with the following parameters:
A-B: check_period = C*10
C-E: check_period = 100.
