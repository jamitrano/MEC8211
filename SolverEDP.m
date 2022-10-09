function [result] = SolverEDP(N,finalTime,numberOfTimeIter,diameter,ratioCoeff,reactionConstant,sourceTerm,dirichletCondition,newmannBorderCondition)
%%% Calcul généraux 
stateVector = zeros(N,1); % vecteur d'etat 
dt = finalTime/numberOfTimeIter;
dx = diameter/2/N;
Dx=FirstDerivateSpaceMatrix(N,dx);
Dxx = SecondDerivateSpaceMatrix(N,dx);
sourceTerm = sourceTerm .*ones(N,1);
[rightMemberMatrix] = DifferentialEquationRightMember(Dx,Dxx,N,dx,ratioCoeff,reactionConstant);

%% Loop 
t = 0:dt:finalTime; %s
nbIter = length(t);
result = zeros(length(stateVector),nbIter);
for i = 1:nbIter
stateVector = EulerImplicitSolverStep(Dx,rightMemberMatrix,sourceTerm,dt,stateVector,dirichletCondition,newmannBorderCondition);
result(:,i) = stateVector;
end
end

