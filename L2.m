function [normL2] = L2(analyticalSol,numericalSol,diameter,N)
%L2 Summary of this function goes here
%   Detailed explanation goes here

normL2=(sum((analyticalSol-numericalSol).^2)/N).^0.5;
end

