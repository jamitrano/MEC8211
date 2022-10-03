function [dfdx] = FirstDerivateSpaceMatrix(n,dx)
%% INPUT 
% n          | size of state vector (N,1) describing the system 
% dx          | Spatial step for discretization 

%% OUPUT 
% dfdx       | Derivation Matrix : df/dx = dfdx*f  
e = ones(n,1);
dfdx = 1/dx *spdiags([-e e],0:1,n,n);
end