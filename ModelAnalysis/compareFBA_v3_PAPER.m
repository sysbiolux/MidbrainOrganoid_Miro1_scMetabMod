%% report fluxSum
clc, close all
clearvars -except solverOK
delete fluxSum.txt
diary fluxSum.txt
diary off
toCompAs=1:4
nrPathways=[1,4,5,6]

plotPerPathway=1
% includeOverviewPlot=nrPathways
includeOverviewPlot=[1,4,5,6,10]

load('consistent_model.mat')

plotResOver=[];
plotResOverPathways={};
for counterP=1:numel(nrPathways)
    clear plotCutoff
    switch nrPathways(counterP)
        case 1
            pathway='Glycolysis'
            metList={'glc_D[c]','g6p[c]','f6p[c]','fdp[c]','dhap[c]','g3p[c]','13dpg[c]','3pg[c]','2pg[c]','pep[c]','pyr[c]'}
        case 2
            pathway='PPP'
            metList={'ru5p_D[c]','s7p[c]','e4p[c]'}
        case 3
            pathway='TCA_(cytoplasm only)'
            metList={'cit[c]','icit[c]','akg[c]','succ[c]','fum[c]','oaa[c]'}
        case 4
            pathway='Tricarboxylic acid cycle'
            metList={'cit[m]','icit[m]','akg[m]','succoa[m]','succ[m]','fum[m]','mal_L[m]','oaa[m]'}
        case 5
            pathway='Oxidative phosphorylation'
            metList={'nadh[m]','fadh2[m]','focytC[m]','q10h2[m]','atp[m]'}
        case 6
            pathway='Fatty acid oxidation'
            metList={'accoa[c]','accoa[m]','coa[c]','coa[m]'}
        case 7
            pathway='Fatty acid oxidation'
            temp=find(ismember(consistent_model.subSystems,pathway));
            metList=findMetsFromRxns(consistent_model,consistent_model.rxns(temp))';
            temp=find(contains(metList,'coa'));
            metList=metList(temp)
            pathway='Fatty acid oxidation__all__coa'
            plotCutoff=10
        case 8
            % Fatty acid oxidation: all carnitin metabolites
            pathway='Fatty acid oxidation'
            temp=find(ismember(consistent_model.subSystems,pathway));
            metList=findMetsFromRxns(consistent_model,consistent_model.rxns(temp))';
            temp=find(contains(metList,'crn'));
            metList=metList(temp)
            pathway='Fatty acid oxidation__all__crn'
            plotCutoff=10
        case 9
            pathway='Cholesterol'
            metList={'sql[r]','zymst[r]','chsterol[r]','chsterol[c]'}
            %             removed: 'mev_R[c]','5pmev[c]',
        case 10
            %             pathway='Tyrosine metabolism'
            %             temp=find(ismember(consistent_model.subSystems,pathway));
            %             metList=findMetsFromRxns(consistent_model,consistent_model.rxns(temp))';
            pathway='Dopamine / Tyrosine metabolism'
            metList={'phe_L[c]','tyr_L[c]','34dhphe[c]','dopa[c]','34dhpac[c]','34dhpha[c]','34dhpe[c]','tym[c]'}
            % removed:
    end
    
    resAllMets=nan(numel(toCompAs),5+1,numel(metList));
    for counterA=1:numel(toCompAs)
        toCompA=toCompAs(counterA)
        load(['FBAmodel' num2str(toCompA) '.mat'])
        %             FBAsolution
        data1=FBAsolution.v;
        model1=multicell_model;
        
        subSystems1=cell(size(model1.subSystems));
        for counter=1:size(model1.subSystems)
            if ~isempty(model1.subSystems{counter})
                subSystems1(counter)=model1.subSystems{counter};
            else
                subSystems1(counter)={'NaN'};
            end
        end
        
        %% calculate flux sum per metabolite
        model=model1;
        v=data1;
        temp=repmat(v',size(model.S,1),1);
        fluxes=model.S.*temp;
        fluxSumP=full(sum((fluxes>0).*fluxes,2));
        fluxSumN=full(sum((fluxes<0).*fluxes,2));
        temp=[fluxSumP, fluxSumN];
        
        resAll=[];
        for counter2=1:numel(metList)
            res=nan(1,5);
            for counter=1:5
                temp2=find(ismember(model.mets,[cell2mat(metList(counter2)) '_' num2str(counter)]));
                if ~isempty(temp2)
                    res(counter)=temp(temp2,1);
                end
            end
            %                 res
            resAll=[resAll; res];
            %         counterA
            %         counter2
            resAllMets(counterA,:,counter2)=[res, sum(res)];
        end
        %             resAll
        resAllt=table(resAll,'RowNames',metList)
        
    end
    diary on
    pathway
    disp('Flux Sum of the following metabolites:')
    disp(metList)
    disp('Rows:Model1-4; Columns:cluster1-5 + sum(1-5)')
    for counter=1:numel(metList)
        disp(metList(counter))
        resAllMets(:,:,counter)
    end
    diary off
    
    if plotPerPathway
        %%
        resAllMets
        %         nansum(resAllMets(:,6,:),3)
        plotRes=reshape(resAllMets(:,6,:),size(resAllMets,1),size(resAllMets,3))
        if exist('plotCutoff')
            plotCutoff
            temp=find(sum(plotRes)>10);
            plotRes=plotRes(:,temp);
            metList=metList(temp);
        end
        %%
        f=figure('Position', [400 300 700 500])
        bar(plotRes, 'stacked')
        title(pathway)
        l=legend(metList,'fontweight','bold','fontsize',10,'Location','eastoutside')
        %         l.Direction='reverse';
        ylabel('FluxSum [a.u.]','fontweight','bold','fontsize',12)
        
        set(gca,'xtick',1:4,'xticklabel',{'WT-D30';'WT-D60'; 'PD-D30';'PD-D60'},'fontweight','bold','fontsize',12);
        xtickangle(45)
        
        temp=erase(pathway,' ');
        temp=erase(temp,'/');
        % savefig(f,['figures/' temp '.fig'])
        
        %%
        if ismember(nrPathways(counterP),includeOverviewPlot)
            plotResOver=[plotResOver; nansum(plotRes,2)'];
            plotResOverPathways=[plotResOverPathways; pathway];
        end
    end
end

%% pathway overview
f=figure('Position', [400 300 700 500])
plot(plotResOver','.-','MarkerSize',80,'LineWidth',4)
title('Pathway Overview')
legend(plotResOverPathways,'fontweight','bold','fontsize',10,'Location','northwest')
ylabel('FluxSum [a.u.]','fontweight','bold','fontsize',12)
set(gca,'xtick',1:4,'xticklabel',{'WT-D30';'WT-D60'; 'PD-D30';'PD-D60'},'fontweight','bold','fontsize',12);
xtickangle(45)
% exportgraphics(f,['figures/pathwayOverview.png'])
