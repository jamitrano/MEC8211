function [C] = AnalyticSolution(sourceTerm,radius,initialConcentration,Deff,radiusVector)
%% Solution analytique de l'équation ellipitique lorsque le mode stationnaire est atteint
%% INPUT
% sourceTerm           | Scalaire - Terme source constant (mol/s)
% radius               | Scalaire - Rayon du pilier 
% initialConcentration | Scalaire - Concentration initale en sel (mol/m3)
% Deff                 | Scalaire - Coefficient de diffusion 
% radiusVector         |  (1,N) - Position de chaque noeud sur le rayon du pilier  

%% OUPUT 
% C                   | (N,1) -  Concentration en selle dans le pilier a l'état
%                     | stationnaire suivant les noeuds du maillage 
C = sourceTerm*radius*radius/(4*Deff).*((radiusVector'/radius).^2 - 1) + initialConcentration;
end

