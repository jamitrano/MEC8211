function [dfdx,dfdx2,constantTerm] = AddNewmannBorderCondition(dfdx,dfdx2,newmannCondition)
%Prise en compte des conditions de Newmann en modifant les matrices des dérivées 
%% INPUT 
% dfdx | (N,N) Matrice des dérivées premières 
% dfdx2 | (N,N) Matrice des dérivées secondes 
% newmannCondition | (M,2) avec M<=N , tel que pour toute ligne i de la
%                    | matrice newmannCondition(i,1) soit le numero du noeud modififé et (i,2)
%                    | la nouvelle valeur
%% OUTPUT 
% dfdx | (N,N) Matrice des dérivées premières 
% dfdx2 | (N,N) Matrice des dérivées secondes 
% constantTerm | (N,1) Vecteur des flux imposés par les conditions de
% Newmann

    rowNumber = newmannCondition(:,1);
    dfdx(rowNumber,:) = 0;
    dfdx2(rowNumber,:) = 0;
    constantTerm = zeros(size(dfdx,1),1);
    for i= length(rowNumber)
        constantTerm(rowNumber(i)) = newmannCondition(i,2);
    end
end

