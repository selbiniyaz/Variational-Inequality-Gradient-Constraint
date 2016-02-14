function fun = f(p,t)
global d;
%vector=ones(1,size(t,2));
%fun=d*vector;
N_t=size(t,2);
M_x=sum(reshape(p(1,t(1:3,:)),3,N_t))/3; % x coordinates at the triangle center
M_y=sum(reshape(p(2,t(1:3,:)),3,N_t))/3; % y coordinates at the triangle center
%fun=10.*(M_x>0.5)+5*(M_y(2)<0.3);
% and(abs(M_x)<0.3,abs(M_y)<0.3);
fun=(-10)*(and(abs(M_x)<0.3,abs(M_y)<0.3))+20*(and(abs(M_x)>0.3,abs(M_y)>0.3));
end
