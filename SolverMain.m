%% MAIN FILE
close all
clear
% Choisir le mode de simulation
simulationMode = 3;  % 1 : Solution numérique
% 2 : Solution Analytique
% 3 : Comparaison Numérique/ Analytique

% Choix des parametres
N = 100 ; % Number of  Node
rangeNode = linspace(5,1e3,20); % pour comparaison
diameter = 1; %m

radiusVector = (1:N)/N*diameter/2;
ratioCoeff = 10e-10; %m2/s
reactionConstant = 0;%4e-9 ; % 1/s
sourceTerm = 10e-8 ; %mol/m3/s
initialConcentration = 10 ; % mol/m3
newmannBorderCondition = [1,0];
dirichletCondition = [N,initialConcentration];
finalTime = 5e8 ; %s
numberOfTimeIter = 1e2 ;
convCriteria=0.001;


if  (simulationMode == 2 )
    %% Solution Analytique
    [analyticResult] = AnalyticSolution(sourceTerm,diameter/2,initialConcentration,ratioCoeff,radiusVector);
        plot(radiusVector,analyticResult);
end
    if (simulationMode == 1 )
         result = SolverEDP(N,finalTime,numberOfTimeIter,convCriteria,diameter,ratioCoeff,reactionConstant,sourceTerm,dirichletCondition,newmannBorderCondition);
        %% Display
        t = 0:finalTime/numberOfTimeIter:finalTime; %s
        dx = diameter/2/N;
        [x,y ] = meshgrid(t,(1:N).*dx);
        figure;
        fig =mesh(x,y,result);
        title('Titre');
        xlabel('Temps (s)');
        ylabel('Distance');
        zlabel('Concentration');
        figure
        plot(radiusVector,result(:,end))
    end
    if  (simulationMode == 3 )
        lNode = length(rangeNode);
        resultOverDx = cell(lNode);
        vecteurLInf=zeros(length(rangeNode),1);
        vecteurL1=zeros(length(rangeNode),1);
        vecteurL2=zeros(length(rangeNode),1);

        for i = 1:lNode
            i
            N = floor(rangeNode(i));
            dirichletCondition = [N,initialConcentration]; 
            resultOverTime = SolverEDP(N,finalTime,numberOfTimeIter,convCriteria,diameter,ratioCoeff,reactionConstant,sourceTerm,dirichletCondition,newmannBorderCondition);
            radiusVector = ((1:N)/N*diameter/2)';
            [analyticResult] = AnalyticSolution(sourceTerm,diameter/2,initialConcentration,ratioCoeff,radiusVector);
            resultOverDx{i,1} = resultOverTime(:,end);
            resultOverDx{i,2} = analyticResult;

            %Création des vecteurs L_inf, L1 et L2
            vecteurLInf(i)=L_inf(resultOverDx{i,1},resultOverDx{i,2});
            vecteurL1(i)=L1(resultOverDx{i,1},resultOverDx{i,2},diameter,N);
            vecteurL2(i)=L2(resultOverDx{i,1},resultOverDx{i,2},diameter,N);
        end
        
        figure
        hold on
        loglog (rangeNode,vecteurL1, LineWidth=2)
        loglog (rangeNode,vecteurL2,LineWidth=2)
        loglog (rangeNode,vecteurLInf,LineWidth=2)
        xlabel : 
        legend('L1','L2','Linf')
    end