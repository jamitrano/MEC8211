%% Fichier principal de Simulation par Differences Finies
% Cours : MEC8211
% Auteurs : Julien Amitrano & Paul Grosfils
% Date : Oct 2022
% Revision : En cours 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear
%%%%%%%%%%%%%%%%% PARAMETRES %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choisir le mode de simulation
simulationMode =4 ;  % 1 : Solution numérique
                     % 2 : Solution Analytique
                     % 3 : Comparaison Numérique/ Analytique
                     % 4 : MNP
                     % 5 : MMS
ordre=2;
N = 100 ; % Number of  Node
rangeNode = linspace(5,100,5); % pour comparaison (simulationMode 3)
diameter = 1; %
finalTime = 1e9 ; %s
dt = 1e4; %s
convCriteria=1e-4;%1e-2;
MNPcriterion = 1e-7;
Deff = 1e-10; %m2/s
reactionConstant = 4e-9 ; % 1/s
sourceTerm =0; %mol/m3/s
initialConcentration = 10 ; % mol/m3
% Conditions aux limites 
newmannCondition = [1,0];
dirichletCondition = [N,initialConcentration];
radius = diameter*0.5;

% Solution analytique pour MMS
syms t r 
analyticSolution(r,t)=-10.*(radius.^2-r.^2) .*exp(t./finalTime)+10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fsurf(analyticSolution,[0 radius 0 finalTime])

%% Calculs généraux 
radiusVector = linspace(0,radius,N);
numberOfTimeIter = floor(finalTime/dt) ;


switch simulationMode
    case 1 %% Simulation Numérique
        [result,convergence,stationnary] = SolverEDP(N,finalTime,numberOfTimeIter,convCriteria,diameter,Deff,reactionConstant,sourceTerm,dirichletCondition,newmannCondition,ordre,0);
        if ~convergence
            warning ('Il n''y a pas eu convergence')
        end 
        plot(radiusVector,result(:,end))
        %figure
        %plot(radiusVector,stationnary)
        if min(result(:,2))~= result(1,2)
            disp('panic')
        end

 
    case 2 %% Solution Analytique
        [analyticResult] = AnalyticSolution(sourceTerm,radius,initialConcentration,Deff,radiusVector);
        plot(radiusVector,analyticResult);

    case 3 %% Comparaison Analytique / Numérique 
        lNode = length(rangeNode);
        L1ErrorVersusOrder = zeros(lNode,2);
        L2ErrorVersusOrder = zeros(lNode,2);
        LInfErrorVersusOrder = zeros(lNode,2);
        for ordre=1:2 % Boucle sur les ordre de differentation spatiale

            resultOverDx = cell(lNode);
            L1Error=zeros(length(rangeNode),1);
            L2Error=zeros(length(rangeNode),1);
            LinfError=zeros(length(rangeNode),1);

            for i = 1:lNode % Boucle sur le nombre de noeuds 
                disp(i)
                N = floor(rangeNode(i));
                dirichletCondition = [N,initialConcentration];
                [result,convergence,stationnary] = SolverEDP(N,finalTime,numberOfTimeIter,convCriteria,diameter,Deff,reactionConstant,sourceTerm,dirichletCondition,newmannCondition,ordre,0);
                if ~convergence
                    disp("Attention il n'a pas convergence, augmenter le temps max");
                end

                radiusVector = linspace(0,diameter/2,N);
                [analyticResult] = AnalyticSolution(sourceTerm,radius,initialConcentration,Deff,radiusVector);
                resultOverDx{i,3}= stationnary;
                resultOverDx{i,1} = result(:,end);
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

        % Affichage 
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
            disp(errorSpace);
        else
            % Affichage
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
            plot(timeVector,errorTime)
            title('Erreur par la méthode MNP')
            xlabel('Temps (s)');
            ylabel('Erreur de discrétisation temporelle')
          
        end
    case 5 %% MMS
        [analyticSolution,residuMMS] = ComputeMMS(reactionConstant,Deff,analyticSolution);
        %% Verification en temps
        numberNode = 6 ;%floor(N/2);
        timeVector  = linspace(0,finalTime,finalTime/dt+1)';
        MMSSolution = analyticSolution(radiusVector,timeVector)';
        sourceTerm = residuMMS(radiusVector,timeVector)';
        conditionInitiale = MMSSolution(:,1);

        [result,convergence,~] = SolverEDP(N,finalTime,numberOfTimeIter,convCriteria,diameter,Deff,reactionConstant,sourceTerm,dirichletCondition,newmannCondition,ordre,numberNode,conditionInitiale);
        if ~convergence
            warning ('Il n''y a pas eu convergence')
        end 

        % Application solution analytique sur le maillage        
        solutionCompare = MMSSolution(numberNode,1:length(result))';

        % Comparer resultat au résidu 
        [~,L2ErrorTime,~] = ComputeError(result,solutionCompare,length(result));
        disp(L2ErrorTime)

        %% Verification en espace 
        % Application solution analytique sur le maillage 
        timeMMS = 1e5; %s
        [result,convergence,~] = SolverEDP(N,timeMMS,numberOfTimeIter,convCriteria,diameter,Deff,reactionConstant,sourceTerm,dirichletCondition,newmannCondition,ordre,0);
        if ~convergence
            warning ('Il n''y a pas eu convergence')
        end 

        solutionCompare = MMSSolution(:,1+timeMMS/dt);

        % Comparer le resultat 
        [~,L2ErrorSpace,~] = ComputeError(result(:,2),solutionCompare,length(result));
%       [~,L2ErrorSpace,~] = ComputeError(result(2:end,2)-residu(2:end),MMSSolution(2:end),length(residu))
        
        disp(L2ErrorSpace)

end