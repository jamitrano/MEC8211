function [normL1] = L1(analyticalSol,numericalSol,N)
%L1 Summary of this function goes here
%   Detailed explanation goes here
%analyticalSol : Vector containing the solutions at each node following the
%analytical method
%numericalSol : Vector containing the solutions at each node following the
%numerical method

%normL1=sum(dx.*abs(analyticalSol-numericalSol))/diameter;
normL1=sum(abs(analyticalSol-numericalSol))/N;

end
