function [result,convergence] = SolverEDP(N,finalTime,numberOfTimeIter,convCriteria,diameter,ratioCoeff,reactionConstant,sourceTerm,dirichletCondition,newmannBorderCondition)
%%% Calcul généraux 
stateVector = zeros(N,1); % vecteur d'etat 
dt = finalTime/numberOfTimeIter;
dx = diameter/2/N;
Dx=FirstDerivateSpaceMatrix(N,dx);
Dxx = SecondDerivateSpaceMatrix(N,dx)./dx;
sourceTerm = sourceTerm .*ones(N,1);
[rightMemberMatrix] = DifferentialEquationRightMember(Dx,Dxx,N,dx,ratioCoeff,reactionConstant);

%% Loop 
%%Initialisation
t = 0:dt:finalTime; %s
nbIter = length(t);
result = zeros(length(stateVector),nbIter);
i=1;
result(:,i)=AddDirichletBorderCondition(stateVector,rightMemberMatrix,dirichletCondition);
delta(i)=1;

%%
while t(i)<finalTime && delta(i)>convCriteria
    i=i+1;
    stateVector = EulerImplicitSolverStep(Dx,rightMemberMatrix,sourceTerm,dt,stateVector,dirichletCondition,newmannBorderCondition);
    result(:,i) = stateVector;
    delta(i)=max(max(abs(result(:,i-1)./result(:,i)-1)));
    
end
if delta(i)>convCriteria
    convergence=0;
else
    convergence=1;
end
result=result(:,1:i);
end
