function [dfdx] = FirstDerivateSpaceMatrix2ndOrder(n,dx)
%% INPUT 
% n          | size of state vector (N,1) describing the system 
% dx          | Spatial step for discretization 

%% OUPUT 
% dfdx       | Derivation Matrix : df/dx = dfdx*f  
e = ones(n,1);
z=zeros (n,1);
dfdx = 1/(2*dx) *spdiags([-e z e],-1:1,n,n);
end