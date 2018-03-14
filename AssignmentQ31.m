%% Assignment 1 %%
% Question 3
close all
clear all
clc

%% Read the binary files
fid1=fopen('sim1.bin'); % (fluctuating u component)
% fid2=fopen('sim2.bin'); % (fluctuating v component)
% fid3=fopen('sim3.bin'); % (fluctuating w component)
global u
n1=8192;
n2=32;
n3=32;
uraw=fread(fid1,'single');
% vraw=fread(fid2,'single');
% wraw=fread(fid3,'single');
itael=0;
for i=1:n1
 for j=1:n2
 for k=1:n3
 itael=itael+1;
 u(i,j,k)=uraw(itael);
%  v(i,j,k)=vraw(itael);
%  w(i,j,k)=wraw(itael);
 end
 end
end

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
Pr = 10000*10^3 ; % [W] Rated Power
Vcut_in = 4 ; % [m/s] Cut in speed
Vcut_out = 25 ; % [m/s] Cut out speed

omega = 0.673 ; % [rad/s] Constant rotational speed
% THIS ONE CHANGE
V_0 = 8 ; % [m/s] Constant wind speed
rho = 1.225 ; % [kg/m3] air mass density
k_emp = 0.6 ; % empirical value used to calculate W_intermediate

%
delta_t = 0.05 ; % [s]
N =3600 ; % [s]   % so the turbulence structure passes the blades
N_element = length(blade_data) ;


global Theta_pitch Theta_cone Theta_tilt Theta_yaw
Theta_pitch = 0.0 ; % [rad]
Theta_cone = 0 ; % [rad]
Theta_tilt = 0 ; % [rad]
Theta_yaw = 0 ; % [rad]


% grid
global Mx My
n_point=32;
size_grid=200;
dxy=size_grid/n_point;

for i=1:n_point
    Mx(i)=-size_grid/2+i*dxy;
    My(i)=-size_grid/2+i*dxy;
end


%% Initialization %%
V0y = 0 ;
V0z = V_0 ;
Wy = zeros(B,N_element,N) ;
Wz =  zeros(B,N_element,N) ;
a_rem = zeros(B,N_element,N) ; 
pz_TS = zeros(N_element,N) ;
py_TS = zeros(N_element,N);

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
    for b=1:3
        % b
        % loop over each element N_element
        for k=1:(N_element-1)
            % k
            
            if (i < 4095)
                Theta_wing=eval(['Theta_wing',num2str(b)]);
                u_turb = velocity_turbulence(blade_data(k),Theta_wing(i),i);
            else
                u_turb = 0 ;
            end
            
            if k==9 && b==1
            u_turb9(i)=u_turb;
            end
            
            [Vrel_y, Vrel_z] = velocity_compute_turb(u_turb,b, blade_data(k), H, Ls, Wy(b,k,i-1), Wz(b,k,i-1), Theta_wing1(i), Theta_wing2(i), Theta_wing3(i) ) ;
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
            
            pz_TS(k,i) = Lift*cos(phi) + Drag*sin(phi) ; % normal
            py_TS(k,i) = Lift*sin(phi) - Drag*cos(phi) ; % tangential
            if k==9
                pz(i) = Lift*cos(phi) + Drag*sin(phi) ; % normal
                py(i) = Lift*sin(phi) - Drag*cos(phi) ; % tangential
            end
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
          
%             dm(k) = blade_data(k)*py(k) ;
%             dP(k) = omega*dm(k) ;
            
            % W_qs(i) = Wz_qs(i) + Wy_qs(i)
            % tau1 = (1.1/(1-1.3*a))*(R/V_0)
            % tau2 = (0.39-0.26*(r/R)^2)*tau1
            % H = W_qs(i)+k_emp*tau1*(((W_qs(i)-W_qs(i-1)/delta_t)
            % Wint(i) = H + (Wint(i-1)-H)*exp(-delta_t/tau1)
            % W(i) = Wint(i) + (W(i-1)-Wint(i))*exp(-delta_t/tau2)
            
        end
        pz_TS(N_element,i) = 0 ;
        py_TS(N_element,i) = 0 ; 
        
 
%         dm(N_element) = 0 ; 
%         dP(N_element) = 0 ; 
        
        % Sanity check with teacher's results
%         if (i==1000)
%         time(i)
%         figure(1)
%         plot(blade_data(:,1), real(pz)) 
%         figure(2) 
%         plot(blade_data(:,1), real(py)) 
%         end
        
        % thrust computation for each blade
        thrust(b,i) = trapz(blade_data(:,1),real(pz_TS(:,i))) ;
        
        % power computation 
%         Power(b,i) = trapz(blade_data(:,1), real(dP)) ;
        
    end
    Thrust (i)= sum(thrust(:,i)) ;
end

% Power_cum = 3*Power(1) 


%% Plots 
figure(1) 
plot(time, real(pz_TS(9,:)))
xlabel('Time [s]','interpreter','latex',  'FontSize', 12)
ylabel('Normal load to rotor plane at r=65.75m [N]','interpreter','latex',  'FontSize', 12)

figure(2)
plot(time, Thrust)
xlabel('Time [s]', 'interpreter','latex', 'FontSize', 12)
ylabel('Thrust [N]', 'interpreter','latex', 'FontSize', 12)

figure(3)
plot(time,u([1:N],16,16)+8)
hold on 
hline = refline([0 8])
hline.Color = 'r'
legend('Turbulence','Wind speed reference')
set(legend,'FontSize',12)
xlabel('Time [s]', 'interpreter','latex', 'FontSize', 12)
ylabel('Wind speed [m/s]', 'interpreter','latex', 'FontSize', 12)
hold off

%% PSD
Thrust=Thrust(150:end);
pz=pz(150:end);

f_low=1/(N*delta_t);
f_high=0.5/delta_t;
fs=1/delta_t;
f=f_low:0.01:f_high;

figure(4)
[Pyy,f]=pwelch(pz-mean(pz),2048,1024,f,fs);
semilogy(2*pi*f/omega,Pyy)
grid on

xlabel('\omega/\omega_{o}','FontSize',14)
ylabel('PSD [kW^{2}/Hz]','FontSize',14)
% title('Power spectrum density of the normal load at r=65.75 m') ;
set(gca,'FontSize',14) ;

figure(5)
[Pyy,f]=pwelch(Thrust-mean(Thrust),256,128,f,fs);
semilogy(2*pi*f/omega,Pyy)
grid on

xlabel('\omega/\omega_{o}','FontSize',14)
ylabel('PSD [kW^{2}/Hz]','FontSize',14)
% title('Power spectrum density of total thrust')
set(gca,'FontSize',14)


