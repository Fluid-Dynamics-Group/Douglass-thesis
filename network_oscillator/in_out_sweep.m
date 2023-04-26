clear; close all; clc;

%% Initial parameters
% ------------------------------------------------------------------
dt = 0.0025; 
nosc_struct = 4; % number of structural modes

%% Load and clean data
% ------------------------------------------------------------------

% Load backrun data
load('Re_100_kb_0625_50000_75000/structure_inner.dat');

% Load perturbed data
load ('POD_results/amplitude.mat');

%% Mode 1
% ------------------------------------------------------------------

% Seclect mode of interest
mode = kb_0625_mode_1;

% select base mode
mode_amplitude_base = mode{1};

% select which perts get selected for train/test
% test_set = [3 6];
train_set = [2 3 4 5 6 7 8 9];

z_amp_1 = [];
z_dot_amp_1 = [];
z_amp_1_test = [];
z_dot_amp_1_test = [];

for int = train_set %training

    % select amplitude pertibation [2 - 9] (1 is unperturbed)
    mode_amplitude_pert = kb_0625_mode_2{int};
    
    % Seperate fluid and structure modes for baseline and perturbed
    mode_amplitude_base_fluid = mode_amplitude_base(:,nosc_struct+1:end);
    mode_amplitude_base_structure = mode_amplitude_base(:,1:nosc_struct);
    mode_amplitude_pert_fluid = mode_amplitude_pert(:,nosc_struct+1:end);
    mode_amplitude_pert_structure = mode_amplitude_pert(:,1:nosc_struct);
    
    % Append back-run data
    mode_amplitude_base_structure = [structure_inner(:,2:5); mode_amplitude_base_structure];
    mode_amplitude_pert_structure = [structure_inner(:,2:5); mode_amplitude_pert_structure];
    
    %% Pre-allocate Arrays
    % ------------------------------------------------------------------
    
    nt = size(mode_amplitude_base_structure,1);
    
    HTampl_base_struct = zeros(nt,nosc_struct);
    HTphase_base_struct = zeros(nt,nosc_struct);
    HTampl_pert_struct = zeros(nt,nosc_struct);
    HTphase_pert_struct = zeros(nt,nosc_struct);
    
    
    %% Hilbert Transfroms
    % ------------------------------------------------------------------
    
    for i = 1:nosc_struct
        s = mode_amplitude_base_structure(:,i);
        HTampl_base_struct(:,i) = abs(hilbert(s)); 
        HTphase_base_struct(:,i) = angle(hilbert(s));
        s = mode_amplitude_pert_structure(:,i);
        HTampl_pert_struct(:,i) = abs(hilbert(s)); 
        HTphase_pert_struct(:,i) = angle(hilbert(s));
    end
    
    % Remove back-run data from hilbert transformed data
    HTampl_base_struct = HTampl_base_struct(25001:end,:);
    HTphase_base_struct = HTphase_base_struct(25001:end,:);
    
    HTampl_pert_struct = HTampl_pert_struct(25001:end,:);
    HTphase_pert_struct = HTphase_pert_struct(25001:end,:);
    
    % Moving average
    struct_osc_pert = movmean(HTampl_pert_struct.*exp(1i*HTphase_pert_struct)./(HTampl_base_struct.*exp(1i*HTphase_base_struct)),3000);
    
    % truncate first and last 2000 time steps
    struct_osc_pert = struct_osc_pert(2000:end-2000,1:3);
    
    % Create array for both base and perturbed cases of fluid mode-pairs
    fluid_osc_base_temp = [mode_amplitude_base_fluid(:,1) + 1i*mode_amplitude_base_fluid(:,2) ...
                           mode_amplitude_base_fluid(:,3) + 1i*mode_amplitude_base_fluid(:,4) ...
                           mode_amplitude_base_fluid(:,5) + 1i*mode_amplitude_base_fluid(:,6) ...
                           mode_amplitude_base_fluid(:,7) + 1i*mode_amplitude_base_fluid(:,8)];
    
    fluid_osc_pert_temp = [mode_amplitude_pert_fluid(:,1) + 1i*mode_amplitude_pert_fluid(:,2) ...
                           mode_amplitude_pert_fluid(:,3) + 1i*mode_amplitude_pert_fluid(:,4) ...
                           mode_amplitude_pert_fluid(:,5) + 1i*mode_amplitude_pert_fluid(:,6) ...
                           mode_amplitude_pert_fluid(:,7) + 1i*mode_amplitude_pert_fluid(:,8)];
    
    % Moving average
    fluid_osc_pert = movmean(fluid_osc_pert_temp./fluid_osc_base_temp,3000);
    
    % truncate first and last 2000 time steps
    fluid_osc_pert = fluid_osc_pert(2000:end-2000,1:4);
    
    % form z and calculate z_dot (derivative of z)
    z     = [struct_osc_pert fluid_osc_pert];
    z_dot = diff(z)/dt;
    z     = z(2:end,:); % reshape to match z_dot
    
    z_dot_amp_1 = [z_dot_amp_1; z_dot];
    z_amp_1 = [z_amp_1; z];

end

%% Mode 2
% ------------------------------------------------------------------

% Seclect mode of interest
mode = kb_0625_mode_2;

% select base mode
mode_amplitude_base = mode{1};

% select which perts get selected for train/test
% test_set = [4 9];
train_set = [2 3 4 5 6 7 8 9];

z_amp_2 = [];
z_dot_amp_2 = [];
z_amp_2_test = [];
z_dot_amp_2_test = [];

for int = train_set 

    % select amplitude pertibation [2 - 9] (1 is unperturbed)
    mode_amplitude_pert = kb_0625_mode_2{int};
    
    % Seperate fluid and structure modes for baseline and perturbed
    mode_amplitude_base_fluid = mode_amplitude_base(:,nosc_struct+1:end);
    mode_amplitude_base_structure = mode_amplitude_base(:,1:nosc_struct);
    mode_amplitude_pert_fluid = mode_amplitude_pert(:,nosc_struct+1:end);
    mode_amplitude_pert_structure = mode_amplitude_pert(:,1:nosc_struct);
    
    % Append back-run data
    mode_amplitude_base_structure = [structure_inner(:,2:5); mode_amplitude_base_structure];
    mode_amplitude_pert_structure = [structure_inner(:,2:5); mode_amplitude_pert_structure];
    
    %% Pre-allocate Arrays
    % ------------------------------------------------------------------
    
    nt = size(mode_amplitude_base_structure,1);
    
    HTampl_base_struct = zeros(nt,nosc_struct);
    HTphase_base_struct = zeros(nt,nosc_struct);
    HTampl_pert_struct = zeros(nt,nosc_struct);
    HTphase_pert_struct = zeros(nt,nosc_struct);
    
    
    %% Hilbert Transfroms
    % ------------------------------------------------------------------
    
    for i = 1:nosc_struct
        s = mode_amplitude_base_structure(:,i);
        HTampl_base_struct(:,i) = abs(hilbert(s)); 
        HTphase_base_struct(:,i) = angle(hilbert(s));
        s = mode_amplitude_pert_structure(:,i);
        HTampl_pert_struct(:,i) = abs(hilbert(s)); 
        HTphase_pert_struct(:,i) = angle(hilbert(s));
    end
    
    % Remove back-run data from hilbert transformed data
    HTampl_base_struct = HTampl_base_struct(25001:end,:);
    HTphase_base_struct = HTphase_base_struct(25001:end,:);
    
    HTampl_pert_struct = HTampl_pert_struct(25001:end,:);
    HTphase_pert_struct = HTphase_pert_struct(25001:end,:);
    
    % Moving average
    struct_osc_pert = movmean(HTampl_pert_struct.*exp(1i*HTphase_pert_struct)./(HTampl_base_struct.*exp(1i*HTphase_base_struct)),3000);
    
    % truncate first and last 2000 time steps
    struct_osc_pert = struct_osc_pert(2000:end-2000,1:3);
    
    % Create array for both base and perturbed cases of fluid mode-pairs
    fluid_osc_base_temp = [mode_amplitude_base_fluid(:,1) + 1i*mode_amplitude_base_fluid(:,2) ...
                           mode_amplitude_base_fluid(:,3) + 1i*mode_amplitude_base_fluid(:,4) ...
                           mode_amplitude_base_fluid(:,5) + 1i*mode_amplitude_base_fluid(:,6) ...
                           mode_amplitude_base_fluid(:,7) + 1i*mode_amplitude_base_fluid(:,8)];
    
    fluid_osc_pert_temp = [mode_amplitude_pert_fluid(:,1) + 1i*mode_amplitude_pert_fluid(:,2) ...
                           mode_amplitude_pert_fluid(:,3) + 1i*mode_amplitude_pert_fluid(:,4) ...
                           mode_amplitude_pert_fluid(:,5) + 1i*mode_amplitude_pert_fluid(:,6) ...
                           mode_amplitude_pert_fluid(:,7) + 1i*mode_amplitude_pert_fluid(:,8)];
    
    % Moving average
    fluid_osc_pert = movmean(fluid_osc_pert_temp./fluid_osc_base_temp,3000);
    
    % truncate first and last 2000 time steps
    fluid_osc_pert = fluid_osc_pert(2000:end-2000,1:4);
    
    % form z and calculate z_dot (derivative of z)
    z     = [struct_osc_pert fluid_osc_pert];
    z_dot = diff(z)/dt;
    z     = z(2:end,:); % reshape to match z_dot
    
    z_dot_amp_2 = [z_dot_amp_2; z_dot];
    z_amp_2 = [z_amp_2; z];

end

in_deg = [];
out_deg = [];

for p_split = 0.01:0.01:1.00

    % Allocate arrays
    z_test = [];
    z_dot_test = [];
    
    z_train_all = [];
    z_dot_train_all = [];

    z_amp_1_test = [];
    z_amp_2_test = [];
    z_dot_amp_1_test = [];
    z_dot_amp_2_test = [];

    % Create test/train split
    [z_amp_1_train, z_amp_1_test] = TestTrainSplit(z_amp_1, p_split);
    [z_amp_2_train, z_amp_2_test] = TestTrainSplit(z_amp_2, p_split);
    [z_dot_amp_1_train, z_dot_amp_1_test] = TestTrainSplit(z_dot_amp_1, p_split);
    [z_dot_amp_2_train, z_dot_amp_2_test] = TestTrainSplit(z_dot_amp_2, p_split);
    
    % Create test arrays
    z_test     = [z_amp_1_test; z_amp_2_test];
    z_dot_test = [z_dot_amp_1_test; z_dot_amp_2_test];
    
    % Create aggrigate arrays
    z_train_all     = [z_amp_1_train; z_amp_2_train];
    z_dot_train_all = [z_dot_amp_1_train; z_dot_amp_2_train];
    
    % Select training data
    % ------------------------------------------------------------------
    z_train     = z_train_all;
    z_dot_train = z_dot_train_all;
    
    % Create Adjacency Matrix
    A = regression_reconstruct(size(z_dot_train ,1), size(z_dot_train ,2), z_train, z_dot_train);

    % Calulate and store in- and out-degree
    in = sum(abs(A),2);
    out = sum(abs(A),1)';

    in_deg = [in_deg, in];
    out_deg = [out_deg, out];

end

for occl = 1:1:7
    figure;
    plot(in_deg(occl,:)./sum(in_deg,1), '-r');hold on;
    plot(out_deg(occl,:)./sum(out_deg,1), '-k');
    legend('in-degree', 'out-degree');
    ylim([0 1]);
    xlim([1 100]);
    hold off;
end

% save('data_full_sam.mat', 'z_dot_amp_1', 'z_dot_amp_2', 'z_amp_1', 'z_amp_2', 'z_dot_test', 'z_test', 'z_train_all', 'z_dot_train_all')

function [dataTraining, dataTesting] = TestTrainSplit(dataA, p)

    % randomly splits "dataA" in to p% training data and (1-p)% testing data
    % p is the proportion of rows to select for training
    
    N = size(dataA,1);           % total number of rows 
    tf = false(N,1);             % create logical index vector
    tf(1:round(p*N)) = true;
    tf = tf(randperm(N));        % randomise order
    
    dataTraining = dataA(tf,:);
    dataTesting  = dataA(~tf,:);

end
