function [L1Error,L2Error,LinfError] = ComputeError(analyticalSol,numericalSol,N)
%% Calcul des  erreurs L1, L2, Linf entre la solution analytique et numérique 

%% INPUT 
% analyticalSol  | (N,1) - Solution analytique à chaque noeud
% numericalSol   | (N,1) - Solution numérique à chaque noeud
% N              | Nombre de noeuds

%% OUPUT 
% L1Error   |Scalaire - Erreur suivant la norme L1 (moyenne pondérée des écarts absolus) 
% L2Error   | Scalaire  Erreur suivant la norme L2 (moyenne pondérée des écarts euclidiens)
% LinfError | Scalaire - Maximum de le l'écart point à point 
L1Error=sum(abs(analyticalSol-numericalSol))/N;
L2Error=(sum((analyticalSol-numericalSol).^2)/N).^0.5;
LinfError = max(abs(numericalSol-analyticalSol));


