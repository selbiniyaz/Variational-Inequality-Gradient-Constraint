function res = L2_grad( vec, triangles, coordinates)
%% 
% 
% this function computes $\|\max(0, |\nabla vec|^2-1)\|^2_{2,\Omega}$
% 
N_t = size(triangles,2);
s = 0; 
for j = 1:N_t
    tri = triangles(1:3,j);
    area=det([1,1,1; coordinates(:,tri)])/2;   
    s=s+max(0,sum((grad(coordinates(:,tri))*vec(tri)).^2)-1)^2*area;
end
res=s;
end

