function [errorSpace,errorTime] = ComputeMNP(N,finalTime,numberOfTimeIter,convCriteria,diameter,Deff,reactionConstant,sourceTerm,dirichletCondition,newmannCondition,ordre,MNPcriterion)
%COMPUTEMNP Summary of this function goes here
%% Calcul en espace et en temps

% Calcul avec maillage fin
[resultOverTime,~,stationnary] = SolverEDP(N,finalTime,numberOfTimeIter,convCriteria,diameter,Deff,reactionConstant,sourceTerm,dirichletCondition,newmannCondition,ordre,N/2);
% fit
radiusVector = linspace(0,diameter/2,N);
timeVector  = linspace(0,finalTime,length(resultOverTime));
polynomSpace=fit(radiusVector',stationnary,'pchipinterp');
polynomTime=fit(timeVector',resultOverTime,'pchipinterp');
% Calcul de la difference
stateVectorSpace = feval(polynomSpace,radiusVector);
stateVectorTime = feval(polynomTime,timeVector);
[~,L2ErrorSpace,~] = ComputeError(stationnary,stateVectorSpace,N);
[~,L2ErrorTime,~] = ComputeError(resultOverTime,stateVectorTime,length(timeVector));
if ~(L2ErrorSpace<MNPcriterion && L2ErrorTime < MNPcriterion)
    error ="L'erreur entre la courbe extrapolée et les données numériques est trop grande";
else
    %Numerical Solution to Nearby Problem
    [d1,d2] = differentiate(polynomSpace,radiusVector);
    [t1 ] = differentiate(polynomTime,timeVector);
    rInv = 1./(linspace(0,N,N)*diameter/(2*N));
    rInv(1)=0;
    errorSpace = Deff*(d2 + rInv'.*d1)  - reactionConstant.*stateVectorSpace -  sourceTerm .*ones(N,1);
    errorTime = t1 - reactionConstant.*stateVectorTime -  sourceTerm ;
end
end

