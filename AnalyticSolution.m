function [C] = AnalyticSolution(sourceTerm,radius,initialConcentration,Deff,radiusVector)
%% Solution analytique de l'équation ellipitique lorsque le mode stationnaire est atteint

C = sourceTerm*radius*radius/(4*Deff).*((radiusVector'/radius).^2 - 1) + initialConcentration;
end

