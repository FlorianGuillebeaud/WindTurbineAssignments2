function V=velocity_turbulence(r,theta,time)
% r element of the blade, theta position of the blade

global u Mx My


ut=zeros(32,32);
for i=1:32
    for j=1:32
        ut(i,j)=u(time,i,j);
    end
end

y=r*sin(theta);
x=r*cos(theta);

V=interp2(Mx,My,ut,x,y);



end
