# Required Programs:

compress.m
draw_one_cell.m
fastbinorv.m
mutation.m
mutrnd_Dunham.m
mybinornd_mex.c
mybinornd_mex.mexa64 (linux)
mybinornd_mex.mexmaci64 (apple)
mybinornd_mex.mexa64 (windows)
removeZeros.m
save_newborns.m
save_winners
save_run_conditions.m
simulateOneWell_NoCost.m
main_simulateManyWells.m

# Readme

## HOW TO RUN

To run with the top-10 strategy, set the parameter
`gc.newborns_per_adult` on line 30 of the program
main_simulateManyWells.m to 10 and run. Might also
work for top-dog strategy.


## DESCRIPTION OF OUTPUT

The program will generate a results folder containing:
	* gc-struct.mat
	* run_conditions.mat
	* gen[i].mat for each generation [i]
	* nb[i].mat for each generation [i]

### Constants

gc-struct.mat and run_conditions.mat contain constants
that define the experimental and biological properties
of the simulation.

### gen[i].mat

Each gen[i].mat file contains the following variables:
	* DT contains the date and time of the run.
	* gen is the generation number.
	* B has information about byproduct. Specifically,
	  B{j}(k) is the amount of byproduct in community j
	  at timepoint k during maturation.
	* Bio_H has information about the biomass of helpers.
	  Specifically, Bio_H{j}(k) is the biomass of helpers
	  in community j at timepoint k.
	* Bio_M{j}(k) is the biomass of manufacturers in
	  community j at timepoint k.
	* P{j}(k) is the amount of product in community j at
	  timepoint k.
	* R{j}(k) is the amount of resource in community j at
	  timepoint k.
	* fp_manu{j}(k) is the fp value of the kth clonal
	  population of manufacturers in community j at the
	  end of maturation.
	* L_help{j} is the per-cell biomass of helper
	  cells in community j at the end of maturation.
	* L_manu{j}(k) is the per-cell biomass of
	  manufacturer cells of the kth clonal population of
	  manufacturers in community j at the end of
	  maturation.
	* n_genos{j} is the total number number of distinct
	  manufacturer clonal populations in community j
	  at the end of maturation.
	* N_help{j} is the number of helper cells in community
	  j at the end of maturation.
	* N_manu{j}(k) is the number of cells in the kth clonal
	  population of manufacturers in community j at the end
	  of maturation.

Each nb[i].mat file contains the following variables:
	* gen is the generation number.
	* newb_fp_manu{j}(k) is the fp value of the kth clonal
	  population of manufacturers in community j at the
	  beginning of maturation.
	* newb_L_help{j} is the per-cell biomass of helper
	  cells in commuinty j at the beginning of maturation.
	* newb_L_manu{j}(k) is the per-cell biomass of
	  manufacturer cells in the kth population of
	  manufacturers in community j at the beginning of
	  maturation.
	* newb_n_genos{j} is the total number number of
	  distinct manufacturer clonal populations in community
	  j at the beginning of maturation.
	* newb_N_help{j} is the number of helper cells in
	  community j at the beginning of maturation.
	* newb_N_manu{j}(k) is the number of cells in the kth
	  clonal population of manufacturers in community j at the
	  beginning of maturation.
	* win_inds is a list of the selected adult communities
	  from the previous round.
