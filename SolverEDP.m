function [result,convergence,stationnary] =SolverEDP(N,finalTime,numberOfTimeIter,convCriteria,diameter,Deff,reactionConstant,sourceTerm,dirichletCondition,newmannBorderCondition,ordre);
%%% Calcul généraux 
stateVector = double(zeros(N,1)); % vecteur d'etat 
dt = finalTime/numberOfTimeIter;
[rightMemberMatrix,rightMemberMatrixStationnary] = ComputeRightMemberMatrix(diameter,N,Deff,reactionConstant,dt,ordre);

sourceTerm = sourceTerm .*ones(N,1);
stationnary = ComputeStationnary(rightMemberMatrixStationnary,sourceTerm,stateVector,dirichletCondition,newmannBorderCondition,ordre);
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
    stateVector = EulerImplicitSolverStep(rightMemberMatrix,sourceTerm,dt,stateVector,dirichletCondition,newmannBorderCondition,ordre);
    result(:,2) = stateVector;
    delta(i)=max(max(abs(result(:,1)./result(:,2)-1)));
    
end
if delta(i)>convCriteria
    convergence=0;
else
    convergence=1;
end
end
