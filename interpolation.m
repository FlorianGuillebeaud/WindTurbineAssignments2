function [Cl_sta, Cd_sta]=interpolation(k, alpha)

tc_1=importdata('FFA-W3-241_ds.txt'); %241
tc_2=importdata('FFA-W3-301_ds.txt'); %301
tc_3=importdata('FFA-W3-360_ds.txt'); %360
tc_4=importdata('FFA-W3-480_ds.txt'); %480
tc_5=importdata('FFA-W3-600_ds.txt'); %600
tc_6=importdata('cylinder_ds.txt'); %100% CILINDER

if k<3
    Cl_sta=interp1(tc_6(:,1),tc_6(:,2),alpha);
    Cd_sta=interp1(tc_6(:,1),tc_6(:,3),alpha);
    % f_star(k)=interp1(tc_6(:,1),tc_6(:,4),alpha(k));
    % Cl_inv(k)=interp1(tc_6(:,1),tc_6(:,5),alpha(k));
    % Cl_fs(k)=interp1(tc_6(:,1),tc_6(:,6),alpha(k));
end
if k>9
    Cl_sta=interp1(tc_1(:,1),tc_1(:,2),alpha);
    Cd_sta=interp1(tc_1(:,1),tc_1(:,3),alpha);
    % f_star(k)=interp1(tc_1(:,1),tc_1(:,4),alpha(k));
    % Cl_inv(k)=interp1(tc_1(:,1),tc_1(:,5),alpha(k));
    % Cl_fs(k)=interp1(tc_1(:,1),tc_1(:,6),alpha(k));
end
if k>6 && k<10
    Cl_sta_interp=[interp1(tc_2(:,1),tc_2(:,2),alpha) interp1(tc_1(:,1),tc_1(:,2),alpha)];
    Cd_sta_interp=[interp1(tc_2(:,1),tc_2(:,3),alpha) interp1(tc_1(:,1),tc_1(:,3),alpha)];
    table=[30.1 24.1];
    if k==9
        Cl_sta=interp1(table, Cl_sta_interp, 24.89);
        Cd_sta=interp1(table, Cd_sta_interp, 24.89);
    end
    if k==8
        Cl_sta=interp1(table, Cl_sta_interp, 25.32);
        Cd_sta=interp1(table, Cd_sta_interp, 25.32);
    end
    if k==7
        Cl_sta=interp1(table, Cl_sta_interp, 27.81);
        Cd_sta=interp1(table, Cd_sta_interp, 27.81);
    end
end
if k==6
    table=[36 30.1];
    Cl_sta_interp=[interp1(tc_3(:,1),tc_3(:,2),alpha) interp1(tc_2(:,1),tc_2(:,2),alpha)];
    Cd_sta_interp=[interp1(tc_3(:,1),tc_3(:,3),alpha) interp1(tc_2(:,1),tc_2(:,3),alpha)];
    
    Cl_sta=interp1(table, Cl_sta_interp, 32.42);
    Cd_sta=interp1(table, Cd_sta_interp, 32.42);
end
if k==5
    table=[48 36];
    Cl_sta_interp=[interp1(tc_4(:,1),tc_4(:,2),alpha) interp1(tc_3(:,1),tc_3(:,2),alpha)];
    Cd_sta_interp=[interp1(tc_4(:,1),tc_4(:,3),alpha) interp1(tc_3(:,1),tc_3(:,3),alpha)];
    
    Cl_sta=interp1(table, Cl_sta_interp, 43.04);
    Cd_sta=interp1(table, Cd_sta_interp, 43.04);
end
if k<5 && k>2
    Cl_sta_interp=[interp1(tc_6(:,1),tc_6(:,2),alpha) interp1(tc_5(:,1),tc_5(:,2),alpha)];
    Cd_sta_interp=[interp1(tc_6(:,1),tc_6(:,3),alpha) interp1(tc_5(:,1),tc_5(:,3),alpha)];
    table=[100 60];
    if k==4
        Cl_sta=interp1(table, Cl_sta_interp, 61.1);
        Cd_sta=interp1(table, Cd_sta_interp, 61.1);
    end
    if k==3
        Cl_sta=interp1(table, Cl_sta_interp, 86.05);
        Cd_sta=interp1(table, Cd_sta_interp, 86.05);
    end
end

end







