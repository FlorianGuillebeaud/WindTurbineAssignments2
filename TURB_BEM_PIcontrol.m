function [Thrust, Power, Maero, omega, Theta_pitch, time,u_turb9]=TURB_BEM_PIcontrol(H, Ls, R, B, omega0, V_0, rho, delta_t, N, N_element, Theta_pitch0, Theta_cone, Theta_tilt, Theta_yaw, Kk, Ki, Kp, Irotor)
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

%EXCEEDS MATRIX DIMENSIONS

% grid
global Mx My
n_point=32;
size_grid=200;
dxy=size_grid/n_point;

for i=1:n_point
    Mx(i)=-size_grid/2+i*dxy;
    My(i)=-size_grid/2+i*dxy;
end

%Initialization 
global W3_100 W3_60 W3_48 W3_36 W3_30 W3_24 blade_data M_G omega_list 
Theta_pitch_i(1)=Theta_pitch0;
Theta_pitch(1) = Theta_pitch0; % [rad]

omega_max=deg2rad(10); %rad/s
omega_ref=1.08; %rad/s
omega(1) = omega0 ;
Theta_max=deg2rad(35);
Theta_min=deg2rad(0); %He said -2 as a possible value in class

V0y = 0 ;
V0z = V_0 ;
Wy = zeros(B,N_element,N) ;
Wz =  zeros(B,N_element,N) ;
a_rem = zeros(B,N_element,N) ; 

p = zeros(B, N) ; 

global a_12 a_21 a_34 
a_12 = [cos(Theta_tilt) sin(Theta_yaw)*sin(Theta_tilt) -cos(Theta_yaw)*sin(Theta_tilt) ;
    0 cos(Theta_yaw) sin(Theta_yaw) ;
    sin(Theta_tilt) -sin(Theta_yaw)*cos(Theta_tilt) cos(Theta_yaw)*cos(Theta_tilt)] ;

a_21 = a_12' ;

a_34 = [cos(Theta_cone) 0 -sin(Theta_cone) ;
    0 1 0 ;
    sin(Theta_cone) 0 cos(Theta_cone)];

a_43 = a_34' ;


time(1) = 0 ;

% Wind position initialization
Theta_wing1(1) = 0 ; % blade 1
Theta_wing2(1) = Theta_wing1(1) + 2*pi/3 ; % blade 2
Theta_wing3(1) = Theta_wing1(1) + 4*pi/3 ; % blade 3

    %% Loop
for i=2:N

    time(i) = time(i-1) + delta_t ;
    Theta_wing1(i) = Theta_wing1(i-1) + omega(i-1)*delta_t ; % blade 1
    Theta_wing2(i) = Theta_wing1(i) + 2*pi/3 ; % blade 2
    Theta_wing3(i) = Theta_wing1(i) + 4*pi/3 ; % blade 3
    
    % loop over each blade B
    for b=1:3
        % b
        % loop over each element N_element
        for k=1:(N_element-1)

           
                Theta_wing=eval(['Theta_wing',num2str(b)]);
                u_turb = velocity_turbulence(blade_data(k),Theta_wing(i),i);

            
             if k==9 && b==1
                  u_turb9(i)=u_turb;
             end
            
            [Vrel_y, Vrel_z] = velocity_compute_turb(u_turb,b, blade_data(k), H, Ls, Wy(b,k,i-1), Wz(b,k,i-1), Theta_wing1(i), Theta_wing2(i), Theta_wing3(i),omega(i-1),V_0, Theta_cone  ) ;
            
            phi = atan(real(-Vrel_z)/real(Vrel_y)) ;
            alpha = radtodeg(phi - (-degtorad(blade_data(k,3)) + Theta_pitch(i-1))) ;
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
            
            pz_TS(k) = Lift*cos(phi) + Drag*sin(phi) ; % normal
            py_TS(k) = Lift*sin(phi) - Drag*cos(phi) ; % tangential

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
                Wz(b,k,i) = - B*Lift*cos(phi)/(4*pi*rho*blade_data(k)*F*(sqrt(V0y^2+(V0z+u_turb+fg*Wz(b,k,i-1))^2))) ;
                Wy(b,k,i) = - B*Lift*sin(phi)/(4*pi*rho*blade_data(k)*F*(sqrt(V0y^2+(V0z+u_turb+fg*Wz(b,k,i-1))^2))) ;
            end
          dm(k) = blade_data(k)*py_TS(k) ;
            dP(k) = omega(i-1)*dm(k) ;
            
        end
        pz_TS(N_element) = 0 ;
        py_TS(N_element) = 0 ; 
        dm(N_element) = 0 ; 
        dP(N_element) = 0 ;

        % thrust computation for each blade
        thrust(b) = trapz(blade_data(:,1),real(pz_TS)) ;
        % power computation 
        Power_b(b) = trapz(blade_data(:,1), real(dP)) ;
        % torque
        Maero_b(b)=trapz(blade_data(:,1), real(dm)) ;
    end
    Thrust (i-1)= sum(thrust) ;
    Power(i-1)=sum(Power_b);
    Maero(i-1)=sum(Maero_b);
    MG(i-1)=interp1(omega_list,M_G,omega(i-1));

    % PI Control 
    GK=1/(1+Theta_pitch(i-1)/Kk);
    Theta_pitch_p(i)=GK*Kp*(omega(i-1)-omega_ref);
    Theta_pitch_i(i)=Theta_pitch_i(i-1)+GK*Ki*(omega(i-1)-omega_ref)*delta_t;
    Theta_pitch_sp=Theta_pitch_p(i)+Theta_pitch_i(i);
    
    if Theta_pitch_sp>(Theta_pitch(i-1)+omega_max*delta_t)
        Theta_pitch(i)=Theta_pitch(i-1)+omega_max*delta_t;
    elseif Theta_pitch_sp<(Theta_pitch(i-1)-omega_max*delta_t)
        Theta_pitch(i)=Theta_pitch(i-1)-omega_max*delta_t;
    else
        Theta_pitch(i)=Theta_pitch_sp;
    end
    if Theta_pitch_sp>=Theta_max
        Theta_pitch(i)=Theta_max;
    elseif Theta_pitch_sp<=Theta_min
        Theta_pitch(i)=Theta_min;
    end
    
    omega(i)=omega(i-1)+(Maero(i-1)-MG(i-1))/Irotor*delta_t;
end
end