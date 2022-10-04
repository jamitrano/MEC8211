%% Test du code 
N = 100 ; % Nn Node 
diameter = 1; %m
stateVector = zeros(N,1); % vecteur d'etat 
ratioCoeff = 1e-10; %m2/an
reactionConstant = 4e-9 ; % 1/an
newmannBorderCondition = [1,0];
dirichletCondition = [N,10];
dx = diameter/2/N;
dt = 1/24; %an
[rightMemberMatrix,constantTerm] = DifferentialEquationRightMember(stateVector,dx,ratioCoeff,reactionConstant,newmannBorderCondition);

%% Loop 
t = 0:dt:1; %an
nbIter = length(t);
result = zeros(length(stateVector),nbIter);
for i = 1:nbIter
stateVector = EulerImplicitSolverStep(rightMemberMatrix,constantTerm,dt,stateVector,dirichletCondition);
result(:,i) = stateVector;
end

%% Display ( a completer ) 
meshGrid = mesh(result)