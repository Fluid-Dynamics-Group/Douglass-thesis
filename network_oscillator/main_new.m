clear;clc;close all;

%% Initial parameters
% ------------------------------------------------------------------
fs    = 14; % font size
lw    = 2; % line width

dt = 0.0025; 
nosc_struct = 4; % number of structural modes

%% Load and clean data
% ------------------------------------------------------------------
load ('POD_results/amplitude.mat');

mode = kb_0625_mode_1;

% select base mode
mode_amplitude_base = mode{1};

% select mode and amplitude pertibation [2 - 9] (1 is unperturbed)
mode_amplitude_pert = mode{9};

% Seperate fluid and structure modes for baseline and perturbed
mode_amplitude_base_fluid = mode_amplitude_base(:,nosc_struct+1:end);
mode_amplitude_base_structure = mode_amplitude_base(:,1:nosc_struct);
mode_amplitude_pert_fluid = mode_amplitude_pert(:,nosc_struct+1:end);
mode_amplitude_pert_structure = mode_amplitude_pert(:,1:nosc_struct);

% Append back-run data
load('Re_100_kb_0625_50000_75000/structure_inner.dat');

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

%% Advance with RK 4-5
% ------------------------------------------------------------------

A = networked_oscillator_model(size(z,1), size(z,2), z, z_dot); % Create occilator network
L = A - diag(sum(A,2)); % Calculate the Laplacian matrix
vf = @(t,x) L*x; %annonymous function
t = 0:dt:24000*dt;
[~,z_model] = ode45(vf, t, z(1,:));


%% Plotting
% ------------------------------------------------------------------

figure;
subplot(221);
plot(abs(z(:,1)),'k-','linewidth',lw);
hold on;
plot(abs(z_model(:,1)),'r-','linewidth',lw);
subplot(222);
plot(abs(z(:,2)),'k-','linewidth',lw);
hold on;
plot(abs(z_model(:,2)),'r-','linewidth',lw);
subplot(223);
plot(abs(z(:,3)),'k-','linewidth',lw);
hold on;
plot(abs(z_model(:,3)),'r-','linewidth',lw);
print('-depsc','structure_dynamics.eps');

figure;
subplot(221);
plot(abs(z(:,4)),'k-','linewidth',lw);
hold on;
plot(abs(z_model(:,4)),'r-','linewidth',lw);
subplot(222);
plot(abs(z(:,5)),'k-','linewidth',lw);
hold on;
plot(abs(z_model(:,5)),'r-','linewidth',lw);
subplot(223);
plot(abs(z(:,6)),'k-','linewidth',lw);
hold on;
plot(abs(z_model(:,6)),'r-','linewidth',lw);
subplot(224);
plot(abs(z(:,7)),'k-','linewidth',lw);
hold on;
plot(abs(z_model(:,7)),'r-','linewidth',lw);
print('-depsc','fluid_dynamics.eps');

adjacency_plotting(abs(A),7,1);

adjacency_plotting(angle(A),7,2);
