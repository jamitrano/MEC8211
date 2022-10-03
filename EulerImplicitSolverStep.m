function [newStepVector] = EulerImplicitSolverStep(rightMember,constantTerm,dt,stateVector,dirichletCondition)
%% Solver df/dt = rightMember*f + constantTerm ; avec une méthode  d'euler implicte 
%% INPUT 
% rightMember  | (N,N) Partie linéaire de l'équation différentielle
% constantTerm | (N,1) Partie affine de l'quation différentielle
% dt           | Pas temporel de résolution 
% stateVector  | (N,1) vecteur d'etat au temps t 
% dirichletCondition | (M,2) avec M<=N , tel que pour toute ligne i de la
%                    | matrice dirichletCondition(i,1) soit le numero du noeud modififé et (i,2)
%                    | la nouvelle valeur
%% OUTPUT 
% newStateVector | (N,1) vecteur d'etat modifié au temps t + dt 

newStepVector = (rightMember.*dt - eye(size(rightMember)))\(constantTerm.*dt - stateVector);
newStepVector = AddDirichletBorderCondition(newStepVector,dirichletCondition);
end

