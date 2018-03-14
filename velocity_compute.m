%% Compute the relative wind velocity in Coordinate system 4 in y and z direction %%

function [Vrel_y, Vrel_z] = velocity_compute(j, r, H, Ls, Wy, Wz, Theta_wing1, Theta_wing2, Theta_wing3, omega )
        
global a_12 a_21 a_34 V_0 Theta_cone
        rt = [H 0 0] ; 
        rs = a_21*[0 0 -Ls]' ;
        
        if j==1 % blade 1  
            
            a_23_1 = [cos(Theta_wing1) sin(Theta_wing1) 0 ; 
                       -sin(Theta_wing1) cos(Theta_wing1) 0 ;
                       0 0 1] ;
            a_14_1 = a_34*a_23_1*a_12 ;
            a_41_1 = a_14_1' ;
            
                % r vector for blade 1
                % rb = a_41_1*[r 0 0]' ;
                % r = rt'+rs+rb ;
                % Velovity 
                    V0_4 = a_14_1*[0 0 V_0]' ;
                    Vrel_y = V0_4(2) + Wy - omega*r*cos(Theta_cone) ;
                    Vrel_z = V0_4(3) + Wz ;
            
        elseif j==2 % blade 2
            a_23_2 = [cos(Theta_wing2) sin(Theta_wing2) 0 ; 
                      -sin(Theta_wing2) cos(Theta_wing2) 0 ;
                      0 0 1] ;
            a_14_2 = a_34*a_23_2*a_12 ;
            a_41_2 = a_14_2' ;            
         
                % r vector for blade 2
                % rb = a_41_2*[r 0 0]' ;
                % r(:,i) = rt'+rs+rb ;
                % Velovity 
                    V0_4 = a_14_2*[0 0 V_0]' ;
                    Vrel_y = V0_4(2) + Wy - omega*r*cos(Theta_cone) ;
                    Vrel_z = V0_4(3) + Wz ;
            
        else % blade 3 
            a_23_3 = [cos(Theta_wing3) sin(Theta_wing3) 0 ; 
                      -sin(Theta_wing3) cos(Theta_wing3) 0 ;
                      0 0 1] ;
            a_14_3 = a_34*a_23_3*a_12 ;
            a_41_3 = a_14_3' ; 
                % r vector for blade 2
                % rb = a_41_3*[r 0 0]' ;
                % r(:,i) = rt'+rs+rb ;
                % Velovity 
                    V0_4 = a_14_3*[0 0 V_0]' ;
                    Vrel_y = V0_4(2) + Wy - omega*r*cos(Theta_cone) ;
                    Vrel_z = V0_4(3) + Wz ;
             
        end
    
end
