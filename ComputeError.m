function  ComputeError(resultOverDx)
%Tracer les erreurs L1, L2, Linf

% L1Error = 1/N*sum(asb(resultOverDx{:,1} - resultOverDx{:,2} ),1);
% L2Error = sqrt(1/N*sum((numeriqueResult - analyticResult ).^2,1));
LinfError = max(abs(numeriqueResult - analyticResult ));

figure
hold on 
loglog(dxVector,L1Error)
loglog(dxVector,L2Error)
loglog(dxVector,LinfError)
legend('erreur L1','erreur L2','erreur L\infnty')
hold off

