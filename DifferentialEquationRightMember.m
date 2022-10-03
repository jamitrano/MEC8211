function [rightMemberMatrix,constantTerm] = DifferentialEquationRightMember(stateVector,dx,ratioCoeff,newmannCondition)
%% Calcul les termes A et B de l'équation générique df/dt = Af + B via le calcul de matrice de dérivation 
%% INPUT 
% stateVector | (N,1) vecteur d'etat au temps t 
% dx          | Pas spatial de discrétisation 
% ratioCoeff  | coefficent d'amplification de l'équation
%             |différentielle
% newmannCondition   | (M,2) avec M<=N , tel que pour toute ligne i de la
%                    | matrice newmannCondition(i,1) soit le numero du noeud modififé et (i,2)
%                    | la nouvelle valeur
%% OUTPUT 
%rightMemberMatrix | (N,N) terme A
% constantTerm     | (N,1) terme B
%% Calcul généraux 
n = length(stateVector);
Dx=FirstDerivateSpaceMatrix(n,dx);
Dxx = SecondDerivateSpaceMatrix(n,dx);
[Dx,Dxx,constantTerm] = AddNewmannBorderCondition(Dx,Dxx,newmannCondition);

%% Calcul spécifique pour la Loi de Fick axisymétrique 
Rinv = diag(1./linspace(0,n-1,n))./dx;
Rinv(1,1) = 0; % Axysymmetry 
rightMemberMatrix = ratioCoeff.*(Dxx + Rinv*Dx);
constantTerm = Rinv.*ratioCoeff*constantTerm;
end

