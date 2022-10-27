function [rightMember,stateVector] = AddNewmannBorderCondition(rightMember,stateVector,newmannCondition,ordre)
%Prise en compte des conditions de Newmann 
%% INPUT 
% stateVector | (N,1) vecteur d'etat au temps t 
% newmannCondition | (M,2) avec M<=N , tel que pour toute ligne i de la
%                    | matrice newmannCondition(i,1) soit le numero du noeud modififé et (i,2)
%                    | la nouvelle valeur du flux en ce noeud
% rightMember        | (N,N) Matrice de l'équation
% ordre              | 1 ou 2 suivant le schema de differentiation choisie
%% OUTPUT 
% stateVector | (N,1) vecteur d'etat modifié au temps t 
% rightMember        | (N,N) Matrice de l'équation modifié pour prendre en
% compte les conditions de Newmann.

    for i= length(newmannCondition(:,1))
        idx = newmannCondition(i,1);
        if ordre==1 % Schema decentré d'ordre 1
            
            rightMember(idx,1:2) = [1 -1];
        end
        if ordre==2 % Schema de Gear, decentré d'ordre 2
            rightMember(idx,1:3) = [-3 4 -1];
        end
        stateVector(idx) = 0;
   
    end
end

