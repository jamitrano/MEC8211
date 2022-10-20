%% MAIN FILE
close all
clear
% Choisir le mode de simulation
simulationMode =1 ;  % 1 : Solution numérique
ordre=1;
% 2 : Solution Analytique
% 3 : Comparaison Numérique/ Analytique

% Choix des parametres
N = 100 ; % Number of  Node
rangeNode = linspace(5,1e2,20); % pour comparaison
diameter = 1; %m

radiusVector = linspace(0,diameter/2,N);
ratioCoeff = 1e-10; %m2/s
reactionConstant = 0;%4e-9 ; % 1/s
sourceTerm = 1e-8 ; %mol/m3/s
initialConcentration = 10 ; % mol/m3
newmannBorderCondition = [1,0];
dirichletCondition = [N,initialConcentration];
finalTime = 5e9 ; %s
numberOfTimeIter = 1e8 ;
dt=finalTime/numberOfTimeIter;
convCriteria=1e-6;


if  (simulationMode == 2 )
    %% Solution Analytique
    [analyticResult] = AnalyticSolution(sourceTerm,diameter/2,initialConcentration,ratioCoeff,radiusVector);
        plot(radiusVector,analyticResult);
end
    if (simulationMode == 1 )
         [result,convergence] = SolverEDP(N,finalTime,numberOfTimeIter,convCriteria,diameter,ratioCoeff,reactionConstant,sourceTerm,dirichletCondition,newmannBorderCondition,ordre);
        %% Display
        t = 0:dt:size(result,2)*dt-1; %s
        dx = diameter/2/N;
        [x,y ] = meshgrid(t,(1:N).*dx);
        figure(3);
        fig =mesh(x,y,result);
        title('Titre');
        xlabel('Temps (s)');
        ylabel('Distance');
        zlabel('Concentration');
        figure
        plot(radiusVector,result(:,end))
    end
    if  (simulationMode == 3 )
        for ordre=1:2
            lNode = length(rangeNode);
            resultOverDx = cell(lNode);
            vecteurLInf=zeros(length(rangeNode),1);
            vecteurL1=zeros(length(rangeNode),1);
            vecteurL2=zeros(length(rangeNode),1);
    
            for i = 1:lNode
                disp(i)
                N = floor(rangeNode(i));
                dirichletCondition = [N,initialConcentration]; 
                [resultOverTime,convergence(i)] = SolverEDP(N,finalTime,numberOfTimeIter,convCriteria,diameter,ratioCoeff,reactionConstant,sourceTerm,dirichletCondition,newmannBorderCondition,ordre);
                radiusVector = ((1:N)/N*diameter/2)';
                [analyticResult] = AnalyticSolution(sourceTerm,diameter/2,initialConcentration,ratioCoeff,radiusVector);
                resultOverDx{i,1} = resultOverTime(:,end);
                resultOverDx{i,2} = analyticResult;
    
                %Création des vecteurs L_inf, L1 et L2
                
                vecteurLInf(i)=L_inf(resultOverDx{i,1},resultOverDx{i,2});
                vecteurL1(i)=L1(resultOverDx{i,1},resultOverDx{i,2},N);
                vecteurL2(i)=L2(resultOverDx{i,1},resultOverDx{i,2},N);
            end
            
            elementSize=0.5./floor(rangeNode);
            loglog (elementSize,vecteurL1,'r', LineWidth=2)
            hold on
            grid on
            loglog (elementSize,vecteurL2,'b',LineWidth=2)
            loglog (elementSize,vecteurLInf,'k',LineWidth=2)
            set(gca, 'XDir','reverse')
            xlabel ('Distance entre les noeuds (m)')
            ylabel('Erreur')
            legend('L1','L2','Linf')
            
            title('Progression de l''erreur en fonction de la distance des noeuds')
        end
    end