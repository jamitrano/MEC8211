function [analyticSolution,sourceTerm] = ComputeMMS(K,Deff)
%COMPUTEMMS Summary of this function goes here

%% Declaration variable symboliques
syms r t C ;
EDP =  diff(C,t) -  Deff .*(diff(diff(C,r),r) + 1/r *diff(C,r)) + K*C;

%% Define u
R = 1/2;
analyticSolution =  exp(t)*sin(pi*r/R);
%% Apply u to EDP
sourceTerm = matlabFunction(compose(EDP,analyticSolution));
end

