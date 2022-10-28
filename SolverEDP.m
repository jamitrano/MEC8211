function [result,convergence,stationnary] =SolverEDP(N,finalTime,numberOfTimeIter,convCriteria,diameter,Deff,reactionConstant,sourceTerm,dirichletCondition,newmannBorderCondition,ordre,saveTime,conditionsInitiales)
%% Fonction auxilliaire permettant de résoudre l'équation aux dérivées partielle 
%% INPUT
% finalTime          | Temps final max s'il n'a pas convergence
% numberOfTimeIter   | Nombre d'itérations temporelle 
% convCriteria       | Critère de convergence pour la résolution numérique
% diameter           | Scalaire - Diamètre du pilier 
% N                  | Nombre de noeuds
% Deff               | Scalaire - Coefficient de diffusion
% reactionConstant   | Scalaire -  Coefficient d'absortion 
% sourceTerm         | Scalaire - Terme source constant (mol/s)
% dirichletCondition | (M,2) avec M<=N , tel que pour toute ligne i de la
%                    | matrice dirichletCondition(i,1) soit le numero du noeud modififé et (i,2)
%                    | la nouvelle valeur
% newmannCondition   | (M,2) avec M<=N , tel que pour toute ligne i de la
%                    | matrice newmannCondition(i,1) soit le numero du noeud modififé et (i,2)
% ordre              | 1 ou 2 suivant le schema de differentiation choisie
% saveTime           | Scalaire - A 0 si on veut un resultat en espace et
%                    | non en temps
%                    | A M <= N pour avoir une resultat en temps au noeud M

%% OUTPUT 
% result      | Si saveTime = 0, renvoie (N,1) resultat de la simulation au
%             | temps final ou à la convergence. Sinon renvoie (T,1) , resultat de la
%             |simulation en fonction du temps au noeud saveTime
% convergence | 1 si le resultat à converger, 0 sinon 
% stationnary | (N,1) résultat stationnaire en espace

%% Calcul généraux 
if ~exist('conditionsInitiales','var')
    stateVector = double(zeros(N,1)); % vecteur d'etat 
else
    stateVector = conditionsInitiales;
end
dt = finalTime/numberOfTimeIter;
[rightMemberMatrix,rightMemberMatrixStationnary] = ComputeRightMemberMatrix(diameter,N,Deff,reactionConstant,dt,ordre);
stationnary=0;
if  isscalar(sourceTerm)
    sourceTerm = sourceTerm .*ones(N,numberOfTimeIter);
    stationnary = ComputeStationnary(rightMemberMatrixStationnary,sourceTerm(:,1),stateVector,dirichletCondition,newmannBorderCondition,ordre);
end

%% Loop 
%%Initialisation
t = 0:dt:finalTime; %s
%nbIter = length(t);
result = zeros(length(stateVector),2);
[result(:,1),~]=AddDirichletBorderCondition(stateVector,rightMemberMatrix,dirichletCondition);
i=1;
delta(i)=1;
if saveTime ~=0
    resultOverTime = zeros(length(t),1);
    resultOverTime(i) = result(saveTime,1);
end
%%
while t(i)<finalTime && delta(i)>convCriteria
    result(:,1) = result(:,2);
    stateVector = EulerImplicitSolverStep(rightMemberMatrix,sourceTerm(:,i),dt,stateVector,dirichletCondition,newmannBorderCondition,ordre);
    result(:,2) = stateVector;
    if saveTime ~=0
        resultOverTime(i+1) = result(saveTime,2);
    end
    i=i+1;
    delta(i)=max(max(abs(result(:,1)./result(:,2)-1)));
    
end
if delta(i)>convCriteria
    convergence=0;
else
    convergence=1;
end

if saveTime ~=0
    result = resultOverTime (1:i);
end
end
