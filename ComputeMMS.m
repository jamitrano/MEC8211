function [analyticSolution,sourceTerm] = ComputeMMS(K,Deff,analyticSolution)
%COMPUTEMMS Summary of this function goes here

%% Declaration variable symboliques
syms r t C ;
EDP =  diff(C,t) -  Deff .*(diff(diff(C,r),r) + 1/r *diff(C,r)) + K*C;

%% Apply u to EDP
sourceTerm = matlabFunction(compose(EDP,analyticSolution));
end

