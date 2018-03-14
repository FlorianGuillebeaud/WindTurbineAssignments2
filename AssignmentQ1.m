%% Assignment 1 %%

close all
clear all
clc

%% PART 1
Vo=11.4;
R = 89.17 ; % [m] Rotor radius
lambda=8;
rho = 1.225 ; % [kg/m3] air mass density
A=pi*R^2;
omega=linspace(0,1.5,100);
Tetha_pitch=0;
P_rated=10.64*10^6;
Cp=P_rated./(0.5.*rho.*Vo.^3.*A);
M_G=0.5*rho*A*R^3*Cp^3*omega.^2./lambda^3;
P=M_G.*omega;
for i=1:length(P)
    if P(i)>P_rated
        P(i)=P_rated
    end
end
M_G=P./omega;

figure(1)
plot(omega,M_G)

figure(2)
plot(omega,P)

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

%% Global parameters %%
global omega V_0 rho k_emp B R
H = 119 ; % hub height (m)
Ls = 7.1 ; % m
R = 89.17 ; % [m] Rotor radius
B = 3 ; % Number of blades
Pr = 10000*10^3 ; % [W] Rated Power
Vcut_in = 4 ; % [m/s] Cut in speed
Vcut_out = 25 ; % [m/s] Cut out speed

omega = 0.673 ; % [rad/s] Constant rotational speed
V_0 = 8 ; % [m/s] Constant wind speed
rho = 1.225 ; % [kg/m3] air mass density
k_emp = 0.6 ; % empirical value used to calculate W_intermediate

%
delta_t = 0.02 ; % [s]
N = 1000 ; % [s]
N_element = length(blade_data) ;


global Theta_pitch Theta_cone Theta_tilt Theta_yaw
Theta_pitch = 0.0 ; % [rad]
Theta_cone = 0 ; % [rad]
Theta_tilt = 0 ; % [rad]
Theta_yaw = 0 ; % [rad]

%% Initialization %%
V0y = 0 ;
V0z = V_0 ;
Wy = zeros(B,N_element,N) ;
Wz =  zeros(B,N_element,N) ;
a_rem = zeros(B,N_element,N) ; 

p = zeros(B, N) ; 
time(1) = 0 ;

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
    
    % loop over each blade B
    for b=1:B
        % b
        % loop over each element N_element
        for k=1:N_element
            % k
            [Vrel_y, Vrel_z] = velocity_compute(b, blade_data(k), H, Ls, Wy(b,k,i-1), Wz(b,k,i-1), Theta_wing1(i), Theta_wing2(i), Theta_wing3(i) ) ;
            %Vrel_y = -omega*blade_data(k)+Wy(b,k,i-1) ;  
            %Vrel_z = V_0+Wz(b,k,i-1);
            
            phi = atan(real(-Vrel_z)/real(Vrel_y)) ;
            alpha = radtodeg(phi - (-degtorad(blade_data(k,3)) + Theta_pitch)) ;
            % alpha

            % first method (doesn't take into account the dynamic stall 
            thick = [100, 60, 48, 36, 30.1, 24.1] ; 
            
            % Cl interpolation 
            cl1 = interp1(W3_100(:,1), W3_100(:,2), alpha) ;
            cl2 = interp1(W3_60(:,1), W3_60(:,2), alpha) ;
            cl3 = interp1(W3_48(:,1), W3_48(:,2), alpha);
            cl4 = interp1(W3_36(:,1), W3_36(:,2), alpha);
            cl5 = interp1(W3_30(:,1), W3_30(:,2), alpha);
            cl6 = interp1(W3_24(:,1), W3_24(:,2), alpha);
            cl_union = [cl1 cl2 cl3 cl4 cl5 cl6] ; 
            Cl = interp1(thick, cl_union, blade_data(k,4)) ;
           
            
            % Cd interpolation 
            cd1 = interp1(W3_100(:,1), W3_100(:,3), alpha) ;
            cd2 = interp1(W3_60(:,1), W3_60(:,3), alpha) ;
            cd3 = interp1(W3_48(:,1), W3_48(:,3), alpha);
            cd4 = interp1(W3_36(:,1), W3_36(:,3), alpha);
            cd5 = interp1(W3_30(:,1), W3_30(:,3), alpha);
            cd6 = interp1(W3_24(:,1), W3_24(:,3), alpha);
            cd_union = [cd1 cd2 cd3 cd4 cd5 cd6] ; 
            Cd = interp1(thick, cd_union, blade_data(k,4)) ;
             
       
            % second method
                % calculate Cl_inv , fstatic, Cl_fs using interpolation
                % calculate tau = 4*c/Vrel
                % update fs(i) = fstatic + (fs(i-1)-fstatic)*exp(-delta_t/tau)
                % calculate Cl(i) = fs(i)* Cl,inv + (1-fs(i))*Cl,fs
            
            Vrel_abs = sqrt(Vrel_y^2+Vrel_z^2) ;
            Lift = 0.5*rho*Vrel_abs^2*Cl*blade_data(k,2) ;
            Drag = 0.5*rho*Vrel_abs^2*Cd*blade_data(k,2) ;
            
            pz(k) = Lift*cos(phi) + Drag*sin(phi) ; % normal
            py(k) = Lift*sin(phi) - Drag*cos(phi) ; % tangential
            
            % without Yaw, a can be calculate as follow : 
            a = abs(Wz(b,k,i-1))/V_0 ;
            a_rem(b,k,i) = a ;
            % with yaw : need to be implemented 
            if a<=1/3
                fg = 1 ;
            else
                fg = 1/4*(5-3*a) ;
            end
            
            % Prand
            f = (B/2)*(R-blade_data(k))/(blade_data(k)*abs(sin(phi)));
            F= 2*acos(exp(-f))/pi;
             
           
            % We add this if statement otherwise the last element if NaN 
            % (F = 0 !! )
            if k==N_element
                Wz(b,k,i) = 0 ; 
                Wy(b,k,i) = 0 ; 
            else   
                Wz(b,k,i) = - B*Lift*cos(phi)/(4*pi*rho*blade_data(k)*F*(sqrt(V0y^2+(V0z+fg*Wz(b,k,i-1))^2))) ;
                Wy(b,k,i) = - B*Lift*sin(phi)/(4*pi*rho*blade_data(k)*F*(sqrt(V0y^2+(V0z+fg*Wz(b,k,i-1))^2))) ;
            end
          
            dm(k) = blade_data(k)*py(k) ;
            dP(k) = omega*dm(k) ;
            
            % W_qs(i) = Wz_qs(i) + Wy_qs(i)
            % tau1 = (1.1/(1-1.3*a))*(R/V_0)
            % tau2 = (0.39-0.26*(r/R)^2)*tau1
            % H = W_qs(i)+k_emp*tau1*(((W_qs(i)-W_qs(i-1)/delta_t)
            % Wint(i) = H + (Wint(i-1)-H)*exp(-delta_t/tau1)
            % W(i) = Wint(i) + (W(i-1)-Wint(i))*exp(-delta_t/tau2)
            
        end
        pz(N_element) = 0 ;
        py(N_element) = 0 ; 
        dm(N_element) = 0 ; 
        dP(N_element) = 0 ; 
        
        % Sanity check with teacher's results
        if (i==1000)
        time(i)
        figure(1)
        plot(blade_data(:,1), real(pz)) 
        figure(2) 
        plot(blade_data(:,1), real(py)) 
        end
        
        % thrust computation for each blade
        thrust(b) = trapz(blade_data(:,1),real(pz)) ;
        
        % power computation 
        Power(b) = trapz(blade_data(:,1), real(dP)) ;
        
    end
end

Power_cum = 3*Power(1) 
Thrust = sum(thrust(:,1)) 

%% Plots 
figure(1) 
plot(blade_data(:,1), real(pz))
xlabel('Element position $[m]$','interpreter','latex',  'FontSize', 12)
ylabel('Normal load [N]','interpreter','latex',  'FontSize', 12)

figure(2)
plot(blade_data(:,1), real(py))
xlabel('Element position [m]', 'interpreter','latex', 'FontSize', 12)
ylabel('Tangential load [N]', 'interpreter','latex', 'FontSize', 12)
ylim([-200 700])
global pt 
eps=0.01;
P = p_compute(blade_data(:,1), blade_data(:,2), blade_data(:,3), eps);

