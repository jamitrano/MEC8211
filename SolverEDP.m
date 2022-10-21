function [result,convergence,stationnary] =SolverEDP(N,finalTime,numberOfTimeIter,convCriteria,diameter,Deff,reactionConstant,sourceTerm,dirichletCondition,newmannBorderCondition,ordre);
%%% Calcul généraux 
stateVector = double(zeros(N,1)); % vecteur d'etat 
dt = finalTime/numberOfTimeIter;
dx = diameter/2/(N-1);
rInv = 1./(linspace(0,N,N)*diameter/(2*N));
rInv(1)=0;
if ordre ==1
    backwardTerm = - Deff/(dx*dx).*ones(1,N);
    centralTerm = 2*Deff/(dx*dx) + Deff/dx.*rInv + 1/dt + reactionConstant;
    forwardTerm = -Deff/(dx*dx) -Deff/dx.*rInv;

    
end

if ordre==2
     backwardTerm = - Deff/(dx*dx).*ones(1,N) + Deff/(2*dx).*rInv;
    centralTerm = (2*Deff/(dx*dx)  + 1/dt + reactionConstant).*ones(1,N);
    forwardTerm = -Deff/(dx*dx) -Deff/(2*dx).*rInv;
end
rightMemberMatrix = diag(backwardTerm(2:N),-1) + diag(centralTerm,0) + diag(forwardTerm(1:N-1),1);

 centralTermStationnary = centralTerm - 1/dt ;
rightMemberMatrixStationnary = diag(backwardTerm(2:N),-1) + diag(centralTermStationnary,0) + diag(forwardTerm(1:N-1),1);
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
%result=result(:,1:i);
end
