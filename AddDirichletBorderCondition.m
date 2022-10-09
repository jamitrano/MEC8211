function [stateVector,rightMember] = AddDirichletBorderCondition(stateVector,rightMember,dirichletCondition)
%% Ajoute les conditions de dirichlet en modifiant le vecteur d'etat
%% INPUT 
% stateVector | (N,1) vecteur d'etat au temps t 
% dirichletCondition | (M,2) avec M<=N , tel que pour toute ligne i de la
%                    | matrice dirichletCondition(i,1) soit le numero du noeud modififé et (i,2)
%                    | la nouvelle valeur
%% OUTPUT 
% stateVector | (N,1) vecteur d'etat modifié au temps t 
    for i= size(dirichletCondition,1)
        idx = dirichletCondition(i,1);
        stateVector(idx) = dirichletCondition(i,2);
        rightMember(idx,:) = 0;
        rightMember(idx,idx) = 1;
       
    end
end

