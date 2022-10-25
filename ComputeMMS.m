function [analyticSolution,sourceTerm] = ComputeMMS(K,Deff)
%COMPUTEMMS Summary of this function goes here

%% Declaration variable symboliques
syms r t R Co;
EDP = @(C) diff(C,t) -  Deff .*(diff(diff(C,r),r) + 1/r *diff(C,r)) + K*C;

%% Define u
analyticSolution = Co *exp(t)*sin(pi*r/R);
%% Apply u to EDP
sourceTerm = EDP(analyticSolution);

end

