%% MAIN FILE
close all
clear
% Choisir le mode de simulation
simulationMode =5 ;  % 1 : Solution numérique
                     % 2 : Solution Analytique
                     % 3 : Comparaison Numérique/ Analytique
                     % 4 : MNP
                     % 5 : MMS
ordre=2;
% Choix des parametres
N = 1000 ; % Number of  Node
rangeNode = linspace(5,100,5); % pour comparaison
diameter = 1; %m

finalTime = 1e5 ; %s
dt = 1e3;
convCriteria=1e-7;
MNPcriterion = 1e-7;
Deff = 1e-10; %m2/s
reactionConstant = 4e-9 ; % 1/s
sourceTerm = 0 ; %mol/m3/s
initialConcentration = 10 ; % mol/m3
%% For calculations
radius = diameter*0.5;
radiusVector = linspace(0,radius,N);
newmannCondition = [1,0];
dirichletCondition = [N,initialConcentration];

numberOfTimeIter = floor(finalTime/dt) ;

switch simulationMode
    case 1 %% Simulation Numérique
        [result,convergence,stationnary] = SolverEDP(N,finalTime,numberOfTimeIter,convCriteria,diameter,Deff,reactionConstant,sourceTerm,dirichletCondition,newmannCondition,ordre);
        plot(radiusVector,result(:,end))
        figure
        plot(radiusVector,stationnary)

    case 2 %% Solution Analytique
        [analyticResult] = AnalyticSolution(sourceTerm,radius,initialConcentration,Deff,radiusVector);
        plot(radiusVector,analyticResult);

    case 3
        lNode = length(rangeNode);
        L1ErrorVersusOrder = zeros(lNode,2);
        L2ErrorVersusOrder = zeros(lNode,2);
        LInfErrorVersusOrder = zeros(lNode,2);
        for ordre=1:2

            resultOverDx = cell(lNode);
            L1Error=zeros(length(rangeNode),1);
            L2Error=zeros(length(rangeNode),1);
            LinfError=zeros(length(rangeNode),1);

            for i = 1:lNode
                disp(i)
                N = floor(rangeNode(i));
                dirichletCondition = [N,initialConcentration];
                [resultOverTime,convergence,stationnary] = SolverEDP(N,finalTime,numberOfTimeIter,convCriteria,diameter,Deff,reactionConstant,sourceTerm,dirichletCondition,newmannCondition,ordre);
                if ~convergence
                    disp("Attention il n'a pas convergence, augmenter le temps max");
                end

                radiusVector = linspace(0,diameter/2,N);
                [analyticResult] = AnalyticSolution(sourceTerm,radius,initialConcentration,Deff,radiusVector);
                resultOverDx{i,3}= stationnary;
                resultOverDx{i,1} = resultOverTime(:,end);
                resultOverDx{i,2} = analyticResult;

                %Création des vecteurs L_inf, L1 et L2
                [L1Error(i),L2Error(i),LinfError(i)] = ComputeError(resultOverDx{i,3},resultOverDx{i,2},N);

            end
            L1ErrorVersusOrder(:,ordre) = L1Error;
            L2ErrorVersusOrder(:,ordre) = L2Error;
            LInfErrorVersusOrder(:,ordre) = LinfError;
        end
        elementSize=radius./floor(rangeNode);
        figure


        loglog(elementSize,L1ErrorVersusOrder(:,1))
        grid on
        hold on
        loglog (elementSize,L2ErrorVersusOrder(:,1))
        loglog (elementSize,LInfErrorVersusOrder(:,1))
        loglog (elementSize,L1ErrorVersusOrder(:,2),"--")
        loglog (elementSize,L2ErrorVersusOrder(:,2),"--")
        loglog (elementSize,LInfErrorVersusOrder(:,2),"--")
        set(gca, 'XDir','reverse')
        xlabel ('Distance entre les noeuds (m)')
        ylabel('Erreur (mol/m3)')
        legend('L1-o1','L2 -o1','Linf-o1','L1-o2','L2 -o2','Linf-o2');

        title('Progression de l''erreur en fonction de la distance des noeuds')
        hold off
    case 4 %% MNP
        [errorSpace,errorTime] = ComputeMNP(N,finalTime,numberOfTimeIter,convCriteria,diameter,Deff,reactionConstant,sourceTerm,dirichletCondition,newmannCondition,ordre,MNPcriterion);
        if  isa(errorSpace,'string')
            disp(error);
        else
            figure
            grid on 
            timeVector  = linspace(0,finalTime,length(errorTime));
            loglog(radiusVector,errorSpace)
            title('Erreur par la méthode MNP')
            xlabel('Distance au centre (m)');
            ylabel('Erreur de discrétisation spatiale')
            set(gca, 'XDir','reverse')
            figure
            grid on 
            loglog(timeVector,errorTime)
            title('Erreur par la méthode MNP')
            xlabel('Temps (s)');
            ylabel('Erreur de discrétisation temporelle')
            set(gca, 'XDir','reverse')
        end
    case 5 %% MMS
        [analyticSolution,sourceTerm] = ComputeMMS(reactionConstant,Deff);
end