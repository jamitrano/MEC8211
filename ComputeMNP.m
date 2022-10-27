function [errorSpace,errorTime] = ComputeMNP(N,finalTime,numberOfTimeIter,convCriteria,diameter,Deff,reactionConstant,sourceTerm,dirichletCondition,newmannCondition,ordre,MNPcriterion)
%% Calcul des solutions par la méthode des Problèmes Proches avec une résolution en espace et en temps
%% INPUT 
% N                  | Nombre de noeuds
% finalTime          | Temps final max s'il n'a pas convergence
% numberOfTimeIter   | Nombre d'itérations temporelle 
% convCriteria       | Critère de convergence pour la résolution numérique
% diameter           | Scalaire - Diamètre du pilier 
% Deff               | Scalaire - Coefficient de diffusion 
% reactionConstant   | Scalaire -  Coefficient d'absortion 
% sourceTerm         | Scalaire - Terme source constant (mol/s)
% dirichletCondition | (M,2) avec M<=N , tel que pour toute ligne i de la
%                    | matrice dirichletCondition(i,1) soit le numero du noeud modififé et (i,2)
%                    | la nouvelle valeur
% newmannCondition   | (M,2) avec M<=N , tel que pour toute ligne i de la
%                    | matrice newmannCondition(i,1) soit le numero du noeud modififé et (i,2)
% ordre              | 1 ou 2 suivant le schema de differentiation choisie
% MNPcriterion       | Critère de vérification de l'acceptalibilité de la
%                    | fonction de fit analytique 

%% OUTPUT 
% errorSpace | (N,1) Erreur  à chaque noeud estimée via la MNP au bout du temps final ou
%            |de la convergence 
% errorTime  | (T,1) Erreur à chaque pas de temps pour le noeud central N/2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calcul avec maillage fin
[resultOverTime,~,stationnary] = SolverEDP(N,finalTime,numberOfTimeIter,convCriteria,diameter,Deff,reactionConstant,sourceTerm,dirichletCondition,newmannCondition,ordre,N/2);
%% Fit
radiusVector = linspace(0,diameter/2,N);
timeVector  = linspace(0,finalTime,length(resultOverTime));
polynomSpace=fit(radiusVector',stationnary,'pchipinterp');
polynomTime=fit(timeVector',resultOverTime,'pchipinterp');
%% Calcul de la difference
stateVectorSpace = feval(polynomSpace,radiusVector);
stateVectorTime = feval(polynomTime,timeVector);
[~,L2ErrorSpace,~] = ComputeError(stationnary,stateVectorSpace,N);
[~,L2ErrorTime,~] = ComputeError(resultOverTime,stateVectorTime,length(timeVector));
if ~(L2ErrorSpace<MNPcriterion && L2ErrorTime < MNPcriterion)
    % Non respect du critère de convergence
    errorSpace ="L'erreur entre la courbe extrapolée et les données numériques est trop grande";
else
    % Numerical Solution to Nearby Problem
    [d1,d2] = differentiate(polynomSpace,radiusVector);
    [t1 ] = differentiate(polynomTime,timeVector);
    rInv = 1./(linspace(0,N,N)*diameter/(2*N));
    rInv(1)=0;
    errorSpace = Deff*(d2 + rInv'.*d1)  - reactionConstant.*stateVectorSpace -  sourceTerm .*ones(N,1);
    errorTime = t1 - reactionConstant.*stateVectorTime -  sourceTerm ;
end
end

