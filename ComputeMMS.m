function [analyticSolution,sourceTerm] = ComputeMMS(K,Deff,analyticSolution)
%% Calcul d'une fonction du terme source correspondant à l'ajustement nécessaire pour faire de la solution analytique choisie 
%% Une solution de l'équation (L(u) = sourceTerm)

%% INPUT 
% K                | Coefficient d'absortion 
% Deff             | Coefficient de diffusion 
% analyticSolution | Fonction symbolique f(r,t) devant être solution de
%                  |l'EDP modifiée
%% OUPUT 
% analyticSolution | Fonction symbolique f(r,t) devant être solution de
%                  |l'EDP modifiée
% sourceTerm       | Fonction symbolique f(r,t) representant le terme
%                  | source pour modifier l'EDP initale


% Define syms 
syms r t C ;
% Define EDP
EDP =  diff(C,t) -  Deff .*(diff(diff(C,r),r) + 1/r *diff(C,r)) + K*C;

% Apply u to EDP
sourceTerm = matlabFunction(compose(EDP,analyticSolution));
analyticSolution = matlabFunction(analyticSolution);
end

