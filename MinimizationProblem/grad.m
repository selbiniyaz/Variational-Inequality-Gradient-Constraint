function gradient = grad(vertices) 
% this function calculates gradients of base functions on the triangle
% given by vertices coordinates (gradients of base functions are constant
% on the triangle)
gradient = ([1,1,1;vertices] \ [0, 0; 1,0; 0,1])';

