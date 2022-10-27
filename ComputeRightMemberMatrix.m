function [rightMemberMatrix,rightMemberMatrixStationnary] = ComputeRightMemberMatrix(diameter,N,Deff,reactionConstant,dt,ordre)
%% Calcul de la matrice principale du probleme ( rightMemberMatrix *C(t+1) = C(t) ) ainsi que de cette matrice à l'état stationnaire
%% INPUT 
% diameter           | Scalaire - Diamètre du pilier 
% N                  | Nombre de noeuds
% Deff               | Scalaire - Coefficient de diffusion
% reactionConstant   | Scalaire -  Coefficient d'absortion 
% dt                 | Scalaire - Pas de temps 
% ordre              | 1 ou 2 suivant le schema de differentiation choisie

%% OUPUT
% rightMemberMatrix            | (N,N) matrice principale du probleme
% rightMemberMatrixStationnary | (N,N) matrice principale du probleme à
%                              | l'état stationnaire

dr = diameter/2/(N-1);
dr2 = dr*dr;
rInv = 1./(linspace(0,N,N)*diameter/(2*N));
rInv(1)=0;
%% Calcul des termes remplissant les tri-diagonales
if ordre ==1 % differention spatiale centrée d'ordre 1 
   
    backwardTerm = - Deff/(dr2).*ones(1,N);
    centralTerm = 2*Deff/(dr2) + Deff/dr.*rInv + 1/dt + reactionConstant;
    forwardTerm = -Deff/(dr2) -Deff/dr.*rInv;

    
end

if ordre==2 % Differention spatiale centrée d'ordre 2 
     backwardTerm = - Deff/(dr2).*ones(1,N) + Deff/(2*dr).*rInv;
    centralTerm = (2*Deff/(dr2)  + 1/dt + reactionConstant).*ones(1,N);
    forwardTerm = -Deff/(dr2) -Deff/(2*dr).*rInv;
end
%% Application des termes aux 3 diagonales 
rightMemberMatrix = diag(backwardTerm(1:N-1),-1) + diag(centralTerm,0) + diag(forwardTerm(2:N),1);
%% Adapatation au cas stationnaire 
 centralTermStationnary = centralTerm - 1/dt ;
rightMemberMatrixStationnary = diag(backwardTerm(2:N),-1) + diag(centralTermStationnary,0) + diag(forwardTerm(1:N-1),1);
end

