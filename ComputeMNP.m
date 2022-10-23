function [error] = ComputeMNP(N,finalTime,numberOfTimeIter,convCriteria,diameter,Deff,reactionConstant,sourceTerm,dirichletCondition,newmannCondition,ordre,MNPcriterion)
%COMPUTEMNP Summary of this function goes here
%   Detailed explanation goes here
% Calcul avec maillage fin
[~,~,stationnary] = SolverEDP(N,finalTime,numberOfTimeIter,convCriteria,diameter,Deff,reactionConstant,sourceTerm,dirichletCondition,newmannCondition,ordre);
% fit
radiusVector = linspace(0,diameter/2,N);
p=fit(radiusVector',stationnary,'pchipinterp');
% Calcul de la difference
stateVector = feval(p,radiusVector);
[~,L2Error,~] = ComputeError(stationnary,stateVector,N);
if ~(L2Error<MNPcriterion)
    error ="L'erreur entre la courbe extrapolée et les données numériques est trop grande";
else
    %Numerical Solution to Nearby Problem
    [d1,d2] = differentiate(p,radiusVector);
    rInv = 1./(linspace(0,N,N)*diameter/(2*N));
    rInv(1)=0;
    error = Deff*(d2 + rInv'.*d1) - sourceTerm .*ones(N,1);
end
end

