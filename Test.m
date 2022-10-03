%% Test du code 
N = 100 ; % Nn Node 
stateVector = zeros(N,1); % vecteur d'etat 
ratioCoeff = 1e-1;
newmannBorderCondition = [1,0];
dirichletCondition = [N,3];
dx = 1e-2;
dt = 1e-3;
[rightMemberMatrix,constantTerm] = DifferentialEquationRightMember(stateVector,dx,ratioCoeff,newmannBorderCondition);

%% Loop 
t = 0:dt:1;
nbIter = length(t);
result = zeros(length(stateVector),nbIter);
for i = 1:nbIter
stateVector = EulerImplicitSolverStep(rightMemberMatrix,constantTerm,dt,stateVector,dirichletCondition);
result(:,i) = stateVector;
end

%% Display ( a completer ) 
meshGrid = mesh(result)