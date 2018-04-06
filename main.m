%% Assignment 2 %%

close all
clear all
clc

%% Read Blade and airfoil Data %%
global blade_data
blade_data = xlsread('Blade_data') ;

global W3_100 W3_60 W3_48 W3_36 W3_30 W3_24
W3_100 = importdata('cylinder_ds.txt'); %100% CILINDER
W3_60  = importdata('FFA-W3-600_ds.txt'); %600
W3_48  = importdata('FFA-W3-480_ds.txt'); %480
W3_36  = importdata('FFA-W3-360_ds.txt'); %360
W3_30  = importdata('FFA-W3-301_ds.txt'); %301
W3_24  = importdata('FFA-W3-241_ds.txt'); %241
N_element = length(blade_data) ;

%% Data about WT
H = 119 ; % hub height (m)
Ls = 7.1 ; % m
R = 89.17 ; % [m] Rotor radius
B = 3 ; % Number of blades
Vcut_in = 4 ; % [m/s] Cut in speed
Vcut_out = 25 ; % [m/s] Cut out speed

Theta_pitch = 0.0 ; % [rad]
Theta_cone = 0 ; % [rad]
Theta_tilt = 0 ; % [rad]
Theta_yaw = 0 ; % [rad]
rho = 1.225 ; % [kg/m3] air mass density

%% time data
delta_t = 0.02 ; % [s]
N = 1000 ; % [s]



%% Question 1
%% Optimum generator characteristic
global M_G omega_list
A=pi*R^2;

lambda=8;


V_0 = 7 ; % [m/s] Constant wind speed
omega = lambda*V_0/R ; % [rad/s] Constant rotational speed

%[~, Power, ~,~]=unsteadyBEM(H, Ls, R, B, omega, V_0, rho, delta_t, N, N_element, Theta_pitch, Theta_cone, Theta_tilt, Theta_yaw);


%Cp7=Power(end)/(0.5*rho*V_0^3*A);
Cp7=0.4316;


omega_list=linspace(0,3,100);
M_G=0.5*rho*A*R^3*Cp7.*omega_list.^2./lambda^3;
P=M_G.*omega_list;
P_rated=10.64*10^6;
for i=1:length(P)
    if P(i)>P_rated
        P(i)=P_rated;
        M_G(i)=P(i)/omega_list(i);
    end
end
%% Plot
figure(1)
plot(omega_list,M_G)
xlabel('\omega','interpreter','latex',  'FontSize', 20)
ylabel('M_G','interpreter','latex',  'FontSize', 20)
title('Generator torque caracteristic')
figure(2)
plot(omega_list,P)
xlabel('\omega','interpreter','latex',  'FontSize', 20)
ylabel('Power','interpreter','latex',  'FontSize', 20)
title('Power caracteristic')
%% Question 2
%% Show that the steady result for a constant wind of 7 m/s (below rated) ends in ?=8 and ?p=0°

Kp=1.5;
Ki=0.64;
Kk=deg2rad(14); 
Irotor=1.6*10^8; %kg.m²

omega0=0.7;
Theta_pitch0=deg2rad(5);

%/!\ every angles have to be in rad, except in blade_data
[thrust, Power, Maero, omega, Theta_pitch, time]=unsteadyBEM_PIcontrol(H, Ls, R, B, omega0, V_0, rho, delta_t, N, N_element, Theta_pitch0, Theta_cone, Theta_tilt, Theta_yaw, Kk, Ki, Kp, Irotor);


%% Plots 
figure; 
plot(time, Theta_pitch)
xlabel('time $[s]$','interpreter','latex',  'FontSize', 20)
ylabel('Theta pitch','interpreter','latex',  'FontSize', 20)
title('V0=7m/s with PI controller')
figure;
plot(time(1:end-1), Maero)
xlabel('time $[s]$','interpreter','latex',  'FontSize', 20)
ylabel('Torque','interpreter','latex',  'FontSize', 20)
title('V0=7m/s with PI controller')
figure;
plot(time, omega)
xlabel('time $[s]$','interpreter','latex',  'FontSize', 20)
ylabel('$\omega$','interpreter','latex',  'FontSize', 20)
title('V0=7m/s with PI controller')
%% Calculate the pitch angle to achieve P_rated and V_0=15m/s
V_0=15;
[thrust, Power, Maero, omega, Theta_pitch, time]=unsteadyBEM_PIcontrol(H, Ls, R, B, omega0, V_0, rho, delta_t, N, N_element, Theta_pitch0, Theta_cone, Theta_tilt, Theta_yaw, Kk, Ki, Kp, Irotor);

%% plots
figure;
plot(time(1:end-1), Power)
xlabel('time $[s]$','interpreter','latex',  'FontSize', 12)
ylabel('Power $[W]$','interpreter','latex',  'FontSize', 12)
title('V0=15m/s with PI controller')
figure; 
plot(time, radtodeg(Theta_pitch))
xlabel('time $[s]$','interpreter','latex',  'FontSize', 12)
ylabel('Theta pitch $[deg]$','interpreter','latex',  'FontSize', 12)
title('V0=15m/s with PI controller')
figure;
plot(time, omega)
xlabel('time $[s]$','interpreter','latex',  'FontSize', 12)
ylabel('$\omega$','interpreter','latex',  'FontSize', 12)
title('V0=15m/s with PI controller')
%% Question 3
V_0=7;
[thrust, Power, Maero, omega, Theta_pitch, time, u_turb9]=TURB_BEM_PIcontrol(H, Ls, R, B, omega0, V_0, rho, delta_t, N, N_element, Theta_pitch0, Theta_cone, Theta_tilt, Theta_yaw, Kk, Ki, Kp, Irotor);

%% plots
figure;
plot(time, V_0+u_turb9)
xlabel('time $[s]$','interpreter','latex',  'FontSize', 12)
ylabel('Wind speed $[W]$','interpreter','latex',  'FontSize', 12)
title('V0=7m/s with PI controller and turbulent wind')
figure;
plot(time(1:end-1), Power)
xlabel('time $[s]$','interpreter','latex',  'FontSize', 12)
ylabel('Power $[W]$','interpreter','latex',  'FontSize', 12)
title('V0=7m/s with PI controller and turbulent wind')
figure; 
plot(time, radtodeg(Theta_pitch))
xlabel('time $[s]$','interpreter','latex',  'FontSize', 12)
ylabel('Theta pitch $[deg]$','interpreter','latex',  'FontSize', 12)
title('V0=7m/s with PI controller and turbulent wind')
figure;
plot(time, omega)
xlabel('time $[s]$','interpreter','latex',  'FontSize', 12)
ylabel('$\omega$','interpreter','latex',  'FontSize', 12)
title('V0=7m/s with PI controller and turbulent wind')
%%
V_0=15;
[thrust, Power, Maero, omega, Theta_pitch, time, u_turb9]=TURB_BEM_PIcontrol(H, Ls, R, B, omega0, V_0, rho, delta_t, N, N_element, Theta_pitch0, Theta_cone, Theta_tilt, Theta_yaw, Kk, Ki, Kp, Irotor);

%% plots
figure;
plot(time, V_0+u_turb9)
xlabel('time $[s]$','interpreter','latex',  'FontSize', 12)
ylabel('Wind speed $[W]$','interpreter','latex',  'FontSize', 12)
title('V0=15m/s with PI controller and turbulent wind')
figure;
plot(time(1:end-1), Power)
xlabel('time $[s]$','interpreter','latex',  'FontSize', 12)
ylabel('Power $[W]$','interpreter','latex',  'FontSize', 12)
title('V0=15m/s with PI controller and turbulent wind')
figure; 
plot(time, radtodeg(Theta_pitch))
xlabel('time $[s]$','interpreter','latex',  'FontSize', 12)
ylabel('Theta pitch $[deg]$','interpreter','latex',  'FontSize', 12)
title('V0=15m/s with PI controller and turbulent wind')
figure;
plot(time, omega)
xlabel('time $[s]$','interpreter','latex',  'FontSize', 12)
ylabel('$\omega$','interpreter','latex',  'FontSize', 12)
title('V0=15m/s with PI controller and turbulent wind')


