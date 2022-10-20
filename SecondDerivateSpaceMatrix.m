function [dfdx2] = SecondDerivateSpaceMatrix(n,dx)
%% INPUT 
% n | size of state vector (N,1) describing the system 
% dx          | Spatial step for discretization 

%% OUPUT 
% dfdx2       | Derivation Matrix for 2nd Order in Space : d2f/dx2 = dfdx2*f  
e = ones(n,1);
dfdx2 = 1/(dx*dx) *spdiags([e -2*e e],-1:1,n,n);
end
