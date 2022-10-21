function [rightMemberMatrix] = DifferentialEquationRightMember(Dx,Dxx,n,dx,ratioCoeff,reactionConstant)
%% Calcul les termes A et B de l'équation générique df/dt = Af + B via le calcul de matrice de dérivation 
%% INPUT 
% stateVector | (N,1) vecteur d'etat au temps t 
% dx          | Pas spatial de discrétisation 
% ratioCoeff  | coefficent de diffusion
% reactionConstant   | coefficient d'attenuation
% newmannCondition   | (M,2) avec M<=N , tel que pour toute ligne i de la
%                    | matrice newmannCondition(i,1) soit le numero du noeud modififé et (i,2)
%                    | la nouvelle valeur
%% OUTPUT 
%rightMemberMatrix | (N,N) terme A
% constantTerm     | (N,1) terme B


%% Calcul spécifique pour la Loi de Fick axisymétrique 
Rinv = diag(1./linspace(0,n,n))./dx;
Rinv(1,1) = 0; % Axysymmetry 
rightMemberMatrix = ratioCoeff.*(Dxx + Rinv*Dx) - reactionConstant.*eye(n);
end

