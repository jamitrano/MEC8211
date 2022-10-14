function [normL_inf] = L_inf(analyticalSol,numericalSol)
%L_INF Summary of this function goes here
%   Detailed explanation goes here
%Prend la différence max entre chaque 
normL_inf = max (abs(numericalSol-analyticalSol));
end

