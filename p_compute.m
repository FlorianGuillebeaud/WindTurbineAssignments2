function [P] = p_compute(r, c, beta, eps)
% input : r, beta and c are lists, eps=0.001


global V_0 Theta_pitch omega rho B R W3_100 W3_60 W3_48 W3_36 W3_30 W3_24 blade_data

lambda = omega*R/V_0;

for k = 1:(length(r)-1)

    sigma = c(k)*B/(2*pi*r(k)); %Fraction of annual area covered by blades

    a = 1;
    ap = 1;
    a_new = 0;
    ap_new = 0;
    count = 0;

    while abs(a-a_new)>eps || count<500 
        count=count+1;
        a = a_new;
        ap = ap_new;
        phi = atan(((1-a)*V_0)/((1+ap)*omega*r(k)));
        alpha= radtodeg(phi - Theta_pitch*pi/180 + deg2rad(beta(k)));
        f = (B/2)*(R-r(k))/(r(k)*sin(phi));
        F= 2*acos(exp(-f))/pi;
        
 
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
            
        Cn = Cl*cos(phi)+Cd*sin(phi);
        Ct= Cl*sin(phi)-Cd*cos(phi);
        Cthrust=((1-a)^2*Cn*sigma)/(sin(phi)^2);

        if a<=1/3
            a_new = 1/((4*F*sin(phi)^2)/(sigma*Cn)+1);
        else
            a_star=Cthrust/(4*F*(1-0.25*(5-3*a)*a));
            a_new=0.1*a_star+0.9*a;
        end
    end
    ap_new = 1/((4*F*sin(phi)*cos(phi))/(sigma*Ct)-1);
    %dM(k) = 0.5*rho*B*V_0*(1-a)*omega*r(k)*(1+ap)*c(k)*Ct*r(k)*r(k)/(sin(phi)*cos(phi));
    %dP(k) = omega*dM(k);
    
    Vrel = ((V_0*(1 - a_new))^2 + (omega*r(k)*(1 + ap_new))^2)^0.5 ;
    pn(k) = 0.5*rho*(Vrel^2)*c(k)*Cn ;
    pt(k) = 0.5*rho*(Vrel^2)*c(k)*Ct ;
    
    %pt=L*sin(phi)-D*cos(phi);
    dM(k) = r(k)*B*pt(k);
    dP(k) = omega*dM(k);
end
dP(18) = 0;
pn(18)=0;
pt(18)=0;
figure(1);
hold on
plot(r,pn,'+')
legend('Unsteady BEM code','Steady BEM code')
figure(2);
hold on
plot(r,pt,'+')
legend('Unsteady BEM code','Steady BEM code')
P = trapz(r,dP);
%Cp = P/(0.5*rho*A*(V_0^3));
end