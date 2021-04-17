clear
resultsfolder = 'results/';
start_cycle = 0;
end_cycle = 2e3;
num_cycles = end_cycle;
cycle_duration = 17;
multiplier = 100;
num_wells = 100;
compressWinner = 1;
digits_fp = 5; 
digits_L = 2;
cw_params = [compressWinner, digits_fp, digits_L];
gc.dt = 0.05;
gc.nsteps = ceil(cycle_duration / gc.dt);
gc.b_Hmax = 0.3;
gc.K_HR = 1/5 * multiplier;
gc.K_MR = 1/3 * multiplier;
gc.K_MB = 1/3 * 100 * multiplier;
gc.b_Mmax = 0.7;
gc.c_RH = 10^-4;
gc.c_RM = 10^-4;
gc.c_BM = 1/3;
gc.r_P = 1;
gc.d_M = 3.5 * 10^-3;
gc.d_H = 1.5 * 10^-3;
gc.R_init = 1 * multiplier;
gc.p_mut = 0.002;
gc.frac_null = 0.5;
gc.sp0 = 0.05;
gc.sn0 = 0.067;
gc.g = 0;
gc.newborns_per_adult = 10;

if mod(cycle_duration,gc.dt) ~= 0
    error('cycle duration must be a multiple of dt')
end
time = 0 : gc.dt : cycle_duration;
rng('shuffle');
seeds = randi(10^8,[num_cycles,1]);


ic_fp_manu = [0.13;zeros(100 - 1,1)];
ic_L_manu = [1;zeros(100 - 1,1)];
ic_N_manu = [60 * multiplier ;zeros(100 - 1,1)];
ic_L_help = 1;
ic_N_help = 40 * multiplier;
ic_n_genos = 1;
N0 = ic_L_manu' * ic_N_manu + ic_L_help' * ic_N_help;

newb_fp_manu = cell(1,num_wells);
newb_L_manu = cell(1,num_wells);
newb_N_manu = cell(1,num_wells);
newb_L_help = cell(1,num_wells);
newb_N_help = cell(1,num_wells);
newb_n_genos = cell(1,num_wells);

for i = 1 : num_wells
    newb_fp_manu{i} = ic_fp_manu;
    newb_L_manu{i} = ic_L_manu;
    newb_N_manu{i} = ic_N_manu;
    newb_L_help{i} = ic_L_help;
    newb_N_help{i} = ic_N_help;
    newb_n_genos{i} = ic_n_genos;
end
len_fp_manu_init = length(ic_fp_manu);
wellPlate_fp_manu = {};
wellPlate_L_manu = {};
wellPlate_N_manu = {};
wellPlate_L_help = {};
wellPlate_N_help = {};
wellPlate_Bio_M = {};
wellPlate_Bio_H = {};
wellPlate_R = {};
wellPlate_B = {};
wellPlate_P = {};
wellPlate_n_genos = {};
if start_cycle > 1
    load([resultsfolder 'nb' num2str(start_cycle) '.mat'])
    load([resultsfolder 'run_conditions.mat'])
    num_cycles = end_cycle;
    rng('shuffle')
    seeds = [seeds(1:start_cycle-1); randi(10^8,[num_cycles-start_cycle+1, 1])];
else
    mkdir(resultsfolder)
    save_newborns(newb_fp_manu,newb_L_manu,newb_N_manu,newb_L_help,newb_N_help,newb_n_genos,0,resultsfolder,[]);
end

save_run_conditions(cycle_duration, num_cycles, multiplier, num_wells,...
    cw_params, gc, time, ic_fp_manu, ic_L_manu, ic_N_manu, ic_L_help,...
    ic_N_help, ic_n_genos, resultsfolder, seeds);

for gen = start_cycle : num_cycles
    rng(seeds(gen));
    parfor i = 1 : num_wells
        [wellPlate_fp_manu{i},wellPlate_L_manu{i},wellPlate_N_manu{i},...
            wellPlate_L_help{i},wellPlate_N_help{i},wellPlate_Bio_M{i},...
            wellPlate_Bio_H{i},wellPlate_R{i},wellPlate_B{i},wellPlate_P{i},...
            wellPlate_n_genos{i},gc_out]...
            = simulateOneWell_NoCost(gc,newb_fp_manu{i},newb_L_manu{i},newb_N_manu{i},...
            newb_L_help{i},newb_N_help{i},newb_n_genos{i});
        if ~isequal(gc_out,gc)
            error('You fool! You must never change your global constants')
        end
    end
%     savedata(wellPlate_fp_manu, wellPlate_L_manu, wellPlate_N_manu,...
%         wellPlate_L_help, wellPlate_N_help, wellPlate_Bio_M,...
%         wellPlate_Bio_H, wellPlate_R, wellPlate_B, wellPlate_P, wellPlate_n_genos, gen, resultsfolder, cw_params);
    P_final = zeros(num_wells,1);
    for i = 1 : num_wells
        P_final(i) = wellPlate_P{i}(end);
    end
    P_sorted = sortrows([P_final,(1:num_wells)'],1);
    
    newb_fp_manu = cell(1,num_wells);
    newb_L_manu = cell(1,num_wells);
    newb_N_manu = cell(1,num_wells);
    newb_L_help = cell(1,num_wells);
    newb_N_help = cell(1,num_wells);
    newb_n_genos = cell(1,num_wells);
    
    num_wells_filled = 0;
    finished_pipetting = false;
    win_inds = [];
    for i = 1 : num_wells
        windex = round(P_sorted(end - i + 1,2));
        win_inds = [win_inds,windex];
        nD = floor((wellPlate_Bio_M{windex}(end) + wellPlate_Bio_H{windex}(end)) / (100 * multiplier));
        target_inds = num_wells_filled + 1 : num_wells_filled + min(nD, gc.newborns_per_adult);
        if target_inds(end) >= num_wells
            target_inds(target_inds > num_wells) = [];
            finished_pipetting = true;
        end
        num_wells_filled = num_wells_filled + min(nD, gc.newborns_per_adult);
        
        [newb_fp_manu,newb_L_manu,newb_N_manu,newb_L_help,newb_N_help,newb_n_genos] = pipette(...
            newb_fp_manu,newb_L_manu,newb_N_manu,newb_L_help,newb_N_help,newb_n_genos,...
            wellPlate_fp_manu{windex},wellPlate_L_manu{windex},wellPlate_N_manu{windex},...
            wellPlate_L_help{windex},wellPlate_N_help{windex}, wellPlate_n_genos{windex},...
            target_inds,nD,cw_params,len_fp_manu_init,N0);
        
        if finished_pipetting
            break
        end
    end
    save_winners(wellPlate_fp_manu, wellPlate_L_manu, wellPlate_N_manu,...
        wellPlate_L_help, wellPlate_N_help, wellPlate_Bio_M,...
        wellPlate_Bio_H, wellPlate_R, wellPlate_B, wellPlate_P, wellPlate_n_genos, gen, resultsfolder, cw_params, win_inds);
    
    save_newborns(newb_fp_manu,newb_L_manu,newb_N_manu,newb_L_help,newb_N_help,newb_n_genos,gen,resultsfolder,win_inds);
end