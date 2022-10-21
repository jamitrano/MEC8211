%% MAIN FILE
close all
clear
% Choisir le mode de simulation
simulationMode =3 ;  % 1 : Solution numérique
% 2 : Solution Analytique
% 3 : Comparaison Numérique/ Analytique
ordre=2;
% Choix des parametres
N = 20 ; % Number of  Node
rangeNode = linspace(5,1e2,10); % pour comparaison
diameter = 1; %m

finalTime = 5e9 ; %s
dt = 1e4;
convCriteria=1e-7;
ratioCoeff = 1e-10; %m2/s
reactionConstant = 0;%4e-9 ; % 1/s
sourceTerm = 1e-8 ; %mol/m3/s
initialConcentration = 10 ; % mol/m3
%% For calculations
radius = diameter*0.5;
radiusVector = linspace(0,radius,N);
newmannBorderCondition = [1,0];
dirichletCondition = [N,initialConcentration];

numberOfTimeIter = floor(finalTime/dt) ;

switch simulationMode
    case 1 %% Simulation Numérique
        [result,convergence,stationnary] = SolverEDP(N,finalTime,numberOfTimeIter,convCriteria,diameter,ratioCoeff,reactionConstant,sourceTerm,dirichletCondition,newmannBorderCondition,ordre);
        plot(radiusVector,result(:,end))
        figure 
        plot(radiusVector,stationnary)

    case 2 %% Solution Analytique
        [analyticResult] = AnalyticSolution(sourceTerm,radius,initialConcentration,ratioCoeff,radiusVector);
        plot(radiusVector,analyticResult);

    case 3
        for ordre=1:2
            lNode = length(rangeNode);
            resultOverDx = cell(lNode);
            L1Error=zeros(length(rangeNode),1);
            L2Error=zeros(length(rangeNode),1);
            LinfError=zeros(length(rangeNode),1);

            for i = 1:lNode
                disp(i)
                N = floor(rangeNode(i));
                dirichletCondition = [N,initialConcentration];
                [resultOverTime,convergence,stationnary] = SolverEDP(N,finalTime,numberOfTimeIter,convCriteria,diameter,ratioCoeff,reactionConstant,sourceTerm,dirichletCondition,newmannBorderCondition,ordre);
                if ~convergence
                    disp(strcat(["Attention, pour ]",N,"[noeuds, il n'a pas convergence, augmenter le temps max"]));
                end
                radiusVector = ((1:N)/N*radius)';
                [analyticResult] = AnalyticSolution(sourceTerm,radius,initialConcentration,ratioCoeff,radiusVector);
                resultOverDx{i,3}= stationnary;
                resultOverDx{i,1} = resultOverTime(:,end);
                resultOverDx{i,2} = analyticResult;

                %Création des vecteurs L_inf, L1 et L2
                [L1Error(i),L2Error(i),LinfError(i)] = ComputeError(resultOverDx{i,3},resultOverDx{i,2},N);
           
            end

            elementSize=radius./floor(rangeNode);
            loglog (elementSize,L1Error,'r', LineWidth=2)
            hold on
            grid on
            loglog (elementSize,L2Error,'b',LineWidth=2)
            loglog (elementSize,LinfError,'k',LineWidth=2)
            set(gca, 'XDir','reverse')
            xlabel ('Distance entre les noeuds (m)')
            ylabel('Erreur')
            legend('L1','L2','Linf')

            title('Progression de l''erreur en fonction de la distance des noeuds')
        end
end