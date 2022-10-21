%% Methode des problèmes proches 

% Solution numérique très raffinée
N = 100 ; % Number of  Node
diameter = 1; %m
radiusVector = (1:N)/N*diameter/2;
ratioCoeff = 10e-10; %m2/s
reactionConstant = 0;%4e-9 ; % 1/s
sourceTerm = 1e-8 ; %mol/m3/s
initialConcentration = 10 ; % mol/m3
newmannBorderCondition = [1,0];
dirichletCondition = [N,initialConcentration];
finalTime = 5e9 ; %s
numberOfTimeIter = 1e4 ;
dt=finalTime/numberOfTimeIter;
convCriteria=0.00001;
ordre = 2;

[result,convergence] = SolverEDP(N,finalTime,numberOfTimeIter,convCriteria,diameter,ratioCoeff,reactionConstant,sourceTerm,dirichletCondition,newmannBorderCondition,ordre);

[m,n]=size(result);
timeVec = 0:dt:finalTime;
timeVec = timeVec(1:n);
% On choisit un noeud 
node = N/2 ;
resOverTime = [timeVec;result(node,:)];
p = polyfit(resOverTime(1,:),resOverTime(2,:),5);
