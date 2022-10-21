function [result,convergence,stationnary] = SolverEDP(N,finalTime,numberOfTimeIter,convCriteria,diameter,ratioCoeff,reactionConstant,sourceTerm,dirichletCondition,newmannBorderCondition,ordre)
%%% Calcul généraux 
stateVector = double(zeros(N,1)); % vecteur d'etat 
dt = finalTime/numberOfTimeIter;
dx = diameter/2/N;

if ordre ==1
    Dx=FirstDerivateSpaceMatrix(N,dx);
end

if ordre==2
    Dx=FirstDerivateSpaceMatrix2ndOrder(N,dx);
end

Dxx = SecondDerivateSpaceMatrix(N,dx)./dx;

sourceTerm = sourceTerm .*ones(N,1);
[rightMemberMatrix] = DifferentialEquationRightMember(Dx,Dxx,N,dx,ratioCoeff,reactionConstant);
stationnary = ComputeStationnary(rightMemberMatrix,sourceTerm,stateVector,dirichletCondition,Dx,dx,newmannBorderCondition,ordre);
%% Loop 
%%Initialisation
t = 0:dt:finalTime; %s
%nbIter = length(t);
result = zeros(length(stateVector),2);
i=1;
delta(i)=1;

%%
while t(i)<finalTime && delta(i)>convCriteria
    i=i+1;
    result(:,1) = result(:,2);
    stateVector = EulerImplicitSolverStep(Dx,dx,rightMemberMatrix,sourceTerm,dt,stateVector,dirichletCondition,newmannBorderCondition,ordre);
    result(:,2) = stateVector;
    delta(i)=max(max(abs(result(:,1)./result(:,2)-1)));
    
end
if delta(i)>convCriteria
    convergence=0;
else
    convergence=1;
end
%result=result(:,1:i);
end
