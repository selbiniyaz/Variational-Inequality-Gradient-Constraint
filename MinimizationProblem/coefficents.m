function c= coefficents(t,Grad_u,gamma)
N_t=size(t,2);
c=zeros(4,N_t);
Grad_u2=Grad_u(1,:).^2+Grad_u(2,:).^2;
Xi=(Grad_u2>=1);
c(1,:)=1+2*gamma*Xi.*(Grad_u2-1+2*Grad_u(1,:).*Grad_u(1,:)); % c_1_1
c(2,:)=4*gamma*Xi.*Grad_u(1,:).*Grad_u(2,:);                % c_2_1
c(3,:)=4*gamma*Xi.*Grad_u(1,:).*Grad_u(2,:);                % c_1_2
c(4,:)=1+2*gamma*Xi.*(Grad_u2-1+2*Grad_u(2,:).*Grad_u(2,:)); % c_2_2
end

