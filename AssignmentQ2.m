% Assignment 1 
% Question 2
close all
clear all
clc

%% Read Blade and airfoil Data %%
blade_data = xlsread('Blade_data') ;

W3_100 = importdata('cylinder_ds.txt'); %100% CILINDER
W3_60  = importdata('FFA-W3-600_ds.txt'); %600
W3_48  = importdata('FFA-W3-480_ds.txt'); %480
W3_36  = importdata('FFA-W3-360_ds.txt'); %360
W3_30  = importdata('FFA-W3-301_ds.txt'); %301
W3_24  = importdata('FFA-W3-241_ds.txt'); %241

%% Global parameters %%
H = 119 ; % hub height (m)
Ls = 7.1 ; % m
R = 89.17 ; % [m] Rotor radius
B = 3 ; % Number of blades
Pr = 10000*10^3 ; % [W] Rated Power
Vcut_in = 4 ; % [m/s] Cut in speed
Vcut_out = 25 ; % [m/s] Cut out speed
global omega V_0 rho k_emp
omega = 0.673 ; % [rad/s] Constant rotational speed
V_0 = 8 ; % [m/s] Constant wind speed
rho = 1.225 ; % [kg/m3] air mass density
k_emp = 0.6 ; % empirical value used to calculate W_intermediate

%
delta_t = 0.02 ; % [s]
N = 10000 ; % [s]
N_element = length(blade_data) ;


global Theta_pitch Theta_cone Theta_tilt Theta_yaw
Theta_pitch = 0 ; % [rad]
Theta_cone = 0 ; % [rad]
Theta_tilt = 0 ; % [rad]
Theta_yaw = 0 ; % [rad]

%% Initialization %%
V0y = 0 ;
V0z = V_0 ;
Wy = zeros(N_element,N) ;
Wz =  zeros(N_element,N) ;
p = zeros(B, N) ; 
time(1) = 0 ;
Wy_qs=zeros(N_element,N);
Wz_qs=zeros(N_element,N);
Wy_int=zeros(N_element,N);
Wz_int=zeros(N_element,N);
Hy=zeros(N_element,N);
Hz=zeros(N_element,N);
fs_time(1)=0;

% Wind position initialization
Theta_wing1(1) = 0 ; % blade 1
Theta_wing2(1) = Theta_wing1(1) + 2*pi/3 ; % blade 2
Theta_wing3(1) = Theta_wing1(1) + 4*pi/3 ; % blade 3

% Rotation matrix %
% considering no roll %
global a_12 a_21 a_34 a_43

a_12 = [cos(Theta_tilt) sin(Theta_yaw)*sin(Theta_tilt) -cos(Theta_yaw)*sin(Theta_tilt) ;
    0 cos(Theta_yaw) sin(Theta_yaw) ;
    sin(Theta_tilt) -sin(Theta_yaw)*cos(Theta_tilt) cos(Theta_yaw)*cos(Theta_tilt)] ;

a_21 = a_12' ;

a_34 = [cos(Theta_cone) 0 -sin(Theta_cone) ;
    0 1 0 ;
    sin(Theta_cone) 0 cos(Theta_cone)];

a_43 = a_34' ;


%% Loop
for i=2:N

    time(i) = time(i-1) + delta_t ;
    Theta_wing1(i) = Theta_wing1(i-1) + omega*delta_t ; % blade 1
    Theta_wing2(i) = Theta_wing1(i) + 2*pi/3 ; % blade 2
    Theta_wing3(i) = Theta_wing1(i) + 4*pi/3 ; % blade 3
    
    %Varying Theta_pitch
    if time(i)>=100 && time(i)<=150
        Theta_pitch =degtorad(2);
    else
        Theta_pitch=0;
    end
    
    % loop over each blade B

    for b=3:3
        % b
        % loop over each element N_element
        for k=1:18
            % k
            [Vrel_y, Vrel_z] = velocity_compute(b, blade_data(k), H, Ls, Wy(k,i-1), Wz(k,i-1), Theta_wing1(i), Theta_wing2(i), Theta_wing3(i) ) ;
            
            phi = atan(real(-Vrel_z)/real(Vrel_y)) ;
            alpha = radtodeg(phi - (-degtorad(blade_data(k,3)) + Theta_pitch)) ;
       
            % second method : we now consider dynamic stall for Cl
            thick = [100, 60, 48, 36, 30.1, 24.1] ;
            
            % Cd interpolation (same as Q1)
            cd1 = interp1(W3_100(:,1), W3_100(:,3), alpha) ;
            cd2 = interp1(W3_60(:,1), W3_60(:,3), alpha) ;
            cd3 = interp1(W3_48(:,1), W3_48(:,3), alpha);
            cd4 = interp1(W3_36(:,1), W3_36(:,3), alpha);
            cd5 = interp1(W3_30(:,1), W3_30(:,3), alpha);
            cd6 = interp1(W3_24(:,1), W3_24(:,3), alpha);
            cd_union = [cd1 cd2 cd3 cd4 cd5 cd6] ; 
            Cd = interp1(thick, cd_union, blade_data(k,4)) ;             
       
            
            cl_inv1 = interp1(W3_100(:,1), W3_100(:,6), alpha) ;
            cl_inv2 = interp1(W3_60(:,1), W3_60(:,6), alpha) ;
            cl_inv3 = interp1(W3_48(:,1), W3_48(:,6), alpha);
            cl_inv4 = interp1(W3_36(:,1), W3_36(:,6), alpha);
            cl_inv5 = interp1(W3_30(:,1), W3_30(:,6), alpha);
            cl_inv6 = interp1(W3_24(:,1), W3_24(:,6), alpha);
            cl_inv_union = [cl_inv1 cl_inv2 cl_inv3 cl_inv4 cl_inv5 cl_inv6] ; 
            Cl_inv = interp1(thick, cl_inv_union, blade_data(k,4)) ; 
            
            fs1 = interp1(W3_100(:,1), W3_100(:,5), alpha) ;
            fs2 = interp1(W3_60(:,1), W3_60(:,5), alpha) ;
            fs3 = interp1(W3_48(:,1), W3_48(:,5), alpha);
            fs4 = interp1(W3_36(:,1), W3_36(:,5), alpha);
            fs5 = interp1(W3_30(:,1), W3_30(:,5), alpha);
            fs6 = interp1(W3_24(:,1), W3_24(:,5), alpha);
            fs_union = [fs1 fs2 fs3 fs4 fs5 fs6] ; 
            fs = interp1(thick, fs_union, blade_data(k,4)) ;    
            
            cl_fs1 = interp1(W3_100(:,1), W3_100(:,7), alpha) ;
            cl_fs2 = interp1(W3_60(:,1), W3_60(:,7), alpha) ;
            cl_fs3 = interp1(W3_48(:,1), W3_48(:,7), alpha);
            cl_fs4 = interp1(W3_36(:,1), W3_36(:,7), alpha);
            cl_fs5 = interp1(W3_30(:,1), W3_30(:,7), alpha);
            cl_fs6 = interp1(W3_24(:,1), W3_24(:,7), alpha);
            cl_fs_union = [cl_fs1 cl_fs2 cl_fs3 cl_fs4 cl_fs5 cl_fs6] ; 
            Cl_fs = interp1(thick, cl_fs_union, blade_data(k,4)) ; 
            
            tau = 4*blade_data(k,2)/sqrt(Vrel_y^2+Vrel_z^2);
            fs_time(i) = fs + (fs_time(i-1)-fs)*exp(-delta_t/tau);
            Cl(i) = fs_time(i)* Cl_inv + (1-fs_time(i))*Cl_fs;
            
            Vrel_abs = sqrt(Vrel_y^2+Vrel_z^2) ;
            Lift = 0.5*rho*Vrel_abs^2*Cl(i)*blade_data(k,2) ;
            Drag = 0.5*rho*Vrel_abs^2*Cd*blade_data(k,2) ;
            pz(k) = Lift*cos(phi) + Drag*sin(phi) ;
            py(k) = Lift*sin(phi) - Drag*cos(phi) ;
            
            % without Yaw, a can be calculate as follow : 
            a = abs(Wz(k,i-1))/V_0;
           
            % with yaw : need to be implemented 
            if a<=1/3
                fg = 1 ;
            else
                fg = 1/4*(5-3*a) ;
            end
            
            
            % Prandl
            f = (B/2)*(R-blade_data(k))/(blade_data(k)*abs(sin(phi)));
            F= 2*acos(exp(-f))/pi;
            
           
            % We add this if statement otherwise the last element if NaN 
            % (F = 0 !! )
            if k==N_element
                Wz_qs(k,i) = 0 ; 
                Wy_qs(k,i) = 0 ; 
            else   
                Wz_qs(k,i) = - B*Lift*cos(phi)/(4*pi*rho*blade_data(k)*F*(sqrt(V0y^2+(V0z+fg*Wz(k,i-1))^2))) ;
                Wy_qs(k,i) = - B*Lift*sin(phi)/(4*pi*rho*blade_data(k)*F*(sqrt(V0y^2+(V0z+fg*Wz(k,i-1))^2))) ;
                
            end
                      
            dm(k) = blade_data(k)*py(k) ;
            dP(k) = omega*dm(k) ;
            
            tau1 = (1.1/(1-1.3*a))*(R/V_0);
            tau2 = (0.39-0.26*(blade_data(k)/R)^2)*tau1;
            
            Hy(k,i) = Wy_qs(k,i)+k_emp*tau1*((Wy_qs(k,i)-Wy_qs(k,i-1))/delta_t);
            Hz(k,i) = Wz_qs(k,i)+k_emp*tau1*((Wz_qs(k,i)-Wz_qs(k,i-1))/delta_t);
            
            Wy_int(k,i) = Hy(k,i) + (Wy_int(k,i-1)-Hy(k,i))*exp(-delta_t/tau1);
            Wz_int(k,i) = Hz(k,i) + (Wz_int(k,i-1)-Hz(k,i))*exp(-delta_t/tau1);
            
            Wy(k,i) = Wy_int(k,i) + (Wy(k,i-1)-Wy_int(k,i))*exp(-delta_t/tau2);
            Wz(k,i) = Wz_int(k,i) + (Wz(k,i-1)-Wz_int(k,i))*exp(-delta_t/tau2);
            
        end
        pz(N_element) = 0 ;
        py(N_element) = 0 ; 
        dm(N_element) = 0 ; 
        dP(N_element) = 0 ; 
        
        % Sanity check with teacher's results
        if i==N
            time(i)
            figure(1)
            plot(blade_data(:,1), real(pz)) 
            figure(2) 
            plot(blade_data(:,1), real(py))
        end
        Thrust(b,i) = trapz(blade_data(:,1),real(pz)) ;
        % power computation 
        Power(b,i) = trapz(blade_data(:,1), real(dP)) ; 
        
    end
%     Power_cum(i) = Power(1,9,i)+Power(2,9,i)+Power(3,9,i) ; 
%     Th_cum(i) = Thrust(1,9,i)+Thrust(2,9,i)+Thrust(3,9,i) ; 
end


%% Question 2 - Explain the new time series for the power and thrust as function of time using 
%e.g. the time history for the induced wind speed at r=65.75 m (k=9)

figure(3)
plot(time, sqrt(Wy(6,:).^2+Wz(6,:).^2))
xlabel('Time $[s]$','interpreter','latex',  'FontSize', 12)
ylabel('Induced wind $[\frac{m}{s}]$','interpreter','latex',  'FontSize', 12)
figure(4)
plot(time,3*Power(3,:)/10^3)
xlabel('Time $[s]$','interpreter','latex',  'FontSize', 12)
ylabel('Power $[kW]$','interpreter','latex',  'FontSize', 12)
hold on 
hline = refline([0 3730])
hline.Color = 'r'
legend('Calculated','DTU 10MW Reference')
hold off
figure(5)
plot(time,3*Thrust(3,:)/10^3)
xlabel('Time $[s]$','interpreter','latex',  'FontSize', 12)
ylabel('Thrust $[kNm]$','interpreter','latex',  'FontSize', 12)
hold on 
hline = refline([0 797])
hline.Color = 'r'
legend('Calculated','DTU 10MW Reference')
hold off
