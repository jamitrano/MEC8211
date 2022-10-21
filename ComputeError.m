function [L1Error,L2Error,LinfError] = ComputeError(analyticalSol,numericalSol,N)
%Tracer les erreurs L1, L2, Linf

L1Error=sum(abs(analyticalSol-numericalSol))/N;
L2Error=(sum((analyticalSol-numericalSol).^2)/N).^0.5;
LinfError = max(abs(numericalSol-analyticalSol));


