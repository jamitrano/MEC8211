function stationnary = ComputeStationnary(rightMemberMatrix,sourceTerm,stateVector,dirichletCondition,newmannCondition,ordre)
%% Calcul de la concentration en fonction des noeuds à l'etat stationnaire 
%% INPUT 
% rightMemberMatrix            | (N,N) matrice principale du probleme
% sourceTerm                   | Scalaire - Terme source constant (mol/s)
% stateVector        | (N,1) - Vecteur d'etat representant le probleme 
% dirichletCondition | (M,2) avec M<=N , tel que pour toute ligne i de la
%                    | matrice dirichletCondition(i,1) soit le numero du noeud modififé et (i,2)
%                    | la nouvelle valeur
% newmannCondition   | (M,2) avec M<=N , tel que pour toute ligne i de la
%                    | matrice newmannCondition(i,1) soit le numero du noeud modififé et (i,2)
% ordre              | 1 ou 2 suivant le schema de differentiation choisie
%% OUPUT 
% stationnary        | (N,1) distribution de la concentration sur les
%                    | noeuds à l'etat stationnaire 


 [rightMember,stateVector] = AddNewmannBorderCondition(rightMemberMatrix,stateVector,newmannCondition,ordre);
 [stateVector,rightMember] = AddDirichletBorderCondition(stateVector,rightMember,dirichletCondition);
 stationnary = rightMember\(-sourceTerm + stateVector);
end

