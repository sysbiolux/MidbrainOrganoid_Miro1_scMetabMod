function [solf,resLactate] = FBA_mediumConc_varO2_v2(modelNr,O2uptake)
%%
load('Claudia_medium_cons_250723.mat')
Results_keep;

%%
sortFlag=1 %1:max-min; 2:max; 3:min
ClusterNames={'c1';'c2'; 'c3'; 'c4';'c5'};
names_col={'Ctrl30';'Ctrl60'; 'PD30'; 'PD60';'GC30'; 'GC60'};
CellNumberPerCluster=[
3193	1692	458	351	1975	960
1196	840	1841	405	1143	348
1083	781	921	412	1121	533
703	573	1292	574	364	228
425	439	461	182	858	719
]

load('simpleRev_recon3.mat')

multicell_model=Results_keep(modelNr).multi_cell_population;

load('Table_medium.mat')
T.medium_rxns_keep(cellfun('isempty',T.medium_rxns_keep))=cellstr('nd');
T=[T; {'Ex_dttp[e]','dttp[e]',1}];
[C,IA,IB]=intersect(multicell_model.rxns,T.medium_rxns_keep);
multicell_model.ub(IA)=T.met_Conc_mM(IB);

multicell_model.ub(find(ismember(multicell_model.rxns,'Ex_lnlc[e]')))=1;

multicell_model.ub(find(ismember(multicell_model.rxns,'Ex_o2s[e]')))=0; %superoxide
multicell_model.ub(find(ismember(multicell_model.rxns,'Ex_o2[e]')))=O2uptake; %21.253*2; %o2, Glc x 2
multicell_model.ub(find(ismember(multicell_model.rxns,'Ex_co2[e]')))=21.253*2; %co2, Glc x 2

for i=1:length(ClusterNames)
    biomass(i) = strcat("biomassmaintenancenoTrTr_",num2str(i));
end 
%for all metabolites that contains [u], R are the reactions that involve
%that metabolite -->
[M,R]=find(multicell_model.S(find(contains(multicell_model.mets,'[u]')),:));
R = unique(R);
list_U_reactions = multicell_model.rxns(R);
for i=1:length(ClusterNames)
    ordre = zeros(1,length(ClusterNames));
    if i==i
        ordre(i) = 1;
    end 
    multicell_model.c(find(contains(multicell_model.rxns,'biomassmaintenance'))) = ordre;
    %Flux Balance Analysis (FBA) is used to find the max optimised value if the
    %model was optimising only for this specific biomass reaction
    FBAsol=optimizeCbModel(multicell_model);
    tenpercent = 0.1*FBAsol.f;
    %set lower bound  
    multicell_model.lb(find(contains(multicell_model.rxns,biomass(i)))) = tenpercent;
end

TotalCells=sum(CellNumberPerCluster);
BiomassFix=CellNumberPerCluster./repmat(TotalCells,size(CellNumberPerCluster,1),1);
sum(BiomassFix);
multicell_model.c(find(contains(multicell_model.rxns,'biomassmaintenance'))) = BiomassFix(:,modelNr);

FBAsolution = optimizeCbModel(multicell_model,'max','zero')

%%
count = 0;
for i=1:length(FBAsolution.v)
    if FBAsolution.v(i)~= 0 
        count = count  + 1;
        indx_non_zero_flux(count) =  i;
    end 
end 

Reactions_nonZero_flux_FBA = multicell_model.rxns(indx_non_zero_flux);

%--> finding the U reactions and searching for the index so that we  can find
%the metabolites present in it with the S matrix 
[M,R]=find(multicell_model.S(find(contains(multicell_model.mets,'[u]')),:));
R = unique(R);
list_U_reactions = multicell_model.rxns(R);
[III,IndxFBA_U_rxns, IndxU] = intersect(Reactions_nonZero_flux_FBA,list_U_reactions);
list_U_reactions_FBAflux = Reactions_nonZero_flux_FBA(IndxFBA_U_rxns);
indx_for_finalmodel = find(ismember(multicell_model.rxns,list_U_reactions_FBAflux));
%--> list of metabolites  that take part in  U reactions and hav e a flux
%(those represent [u]metabolites and other metabolites present in U
%reactions)
[MM,RR]=find(multicell_model.S(:,indx_for_finalmodel));
list_of_metabolites_that_have_flux_and_take_part_in_U_reactions = unique(multicell_model.mets(MM));
%--> list of metabolite that take part in U reactions and that have a flux AND
%those are exclusively [u] metabolites
list_of_U_metabolites_that_have_flux = list_of_metabolites_that_have_flux_and_take_part_in_U_reactions(find(contains(list_of_metabolites_that_have_flux_and_take_part_in_U_reactions,'[u]')));

%%
Flux_non_zero_FBA = FBAsolution.v(indx_non_zero_flux);

PerfectMatrix = array2table(zeros(length(list_of_U_metabolites_that_have_flux),length(list_U_reactions_FBAflux)));
PerfectMatrix.Properties.RowNames =  list_of_U_metabolites_that_have_flux(:,1);
PerfectMatrix.Properties.VariableNames = list_U_reactions_FBAflux(:,1);

for i = 1:length(list_of_U_metabolites_that_have_flux)
    %model index for the specific metabolite (i) 
    Mindx = find(ismember(multicell_model.mets,list_of_U_metabolites_that_have_flux(i)));
    
    %model index of all reaction of the specific metabolite THAT contain a flux
    IndxOfSpecificMet_FinalModel = ismember(multicell_model.mets,list_of_U_metabolites_that_have_flux(i));
    [M,R]=find(multicell_model.S(IndxOfSpecificMet_FinalModel,:));
    [I,IA,IB]= intersect(multicell_model.rxns(R),list_U_reactions_FBAflux);
    Rindx = R(IA); %--> model index of all reactions w/ flux AND w/ that specific met 

    %to note down the flux of the reactions for that specific metabolite
    %decide if the metabolite is consumed or produced depending on which
    %side of the reaction it is (S matrix value) and depending on its flux
    % --> THIS  CODE HAS BEEN CORRECTED THROUGH THE TROUBLE SHOOTING PART
    % (SEE DOWN BELOW) 
    for v=1:length(Rindx)
        %right and positive flux --> produced
        if multicell_model.S(Mindx,Rindx(v)) > 0 && Flux_non_zero_FBA(IB(v)) > 0
            PerfectMatrix(multicell_model.mets(Mindx),multicell_model.rxns(Rindx(v))) = array2table(FBAsolution.v(Rindx(v)));
        %left and positive flux --> consumed
        elseif multicell_model.S(Mindx,Rindx(v)) < 0 && Flux_non_zero_FBA(IB(v)) > 0
            PerfectMatrix(multicell_model.mets(Mindx),multicell_model.rxns(Rindx(v))) = array2table(-FBAsolution.v(Rindx(v)));
        %right and negative flux --> consumed
        elseif multicell_model.S(Mindx,Rindx(v)) > 0 && Flux_non_zero_FBA(IB(v)) < 0
            PerfectMatrix(multicell_model.mets(Mindx),multicell_model.rxns(Rindx(v))) = array2table(FBAsolution.v(Rindx(v)));
        %left and negative flux --> produced
        elseif multicell_model.S(Mindx,Rindx(v)) < 0 && Flux_non_zero_FBA(IB(v)) < 0
            PerfectMatrix(multicell_model.mets(Mindx),multicell_model.rxns(Rindx(v))) = array2table(-FBAsolution.v(Rindx(v)));
        end 
    end 
end 

%Creating a new string for the cluster names 
Cluster_names_with_externalCompartment = ClusterNames;
Cluster_names_with_externalCompartment(length(ClusterNames)+1)= {'External Compartment'};

Cluster_annotation = cell(1,length(Cluster_names_with_externalCompartment));
for e = 1:length(Cluster_names_with_externalCompartment)
    if e < length(Cluster_names_with_externalCompartment)
        Cluster_annotation(e) = cellstr(num2str(e)) ;
    elseif e == length(Cluster_names_with_externalCompartment)
        Cluster_annotation(e)= cellstr("e]");
    end 
end 


XOXO = struct;
LastTwo = cellfun(@(S) S(end-1:end), PerfectMatrix.Properties.VariableNames, 'Uniform', 0);
LastThree = cellfun(@(S) S(end-2:end), PerfectMatrix.Properties.VariableNames, 'Uniform', 0);
count = zeros(1,length(Cluster_names_with_externalCompartment));
 

for g=1:length(PerfectMatrix.Properties.VariableNames)
    for h = 1:length(unique(LastTwo)) %TS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if ismember(LastTwo(g),strcat("_",Cluster_annotation{h})) ==1
           count(h) = count(h) + 1;
           XOXO(h).matrices(:,count(h)) = PerfectMatrix(:,g);
        elseif ismember(LastTwo(g),Cluster_annotation{h}) ==1
            count(h) = count(h) + 1;
           XOXO(h).matrices(:,count(h)) = PerfectMatrix(:,g);
        elseif ismember(LastThree(g),strcat("_",Cluster_annotation{h})) ==1
            count(h) = count(h) + 1;
            XOXO(h).matrices(:,count(h)) = PerfectMatrix(:,g);
        end 
    end 
end 

%Summing up all the columns of each matrix in order to create only one summed up column for each cluster  
Sum_of_XOXO_columns = struct;
for i= 1:numel(XOXO) %width
    for w=1:height(XOXO(i).matrices)
        Sum_of_XOXO_columns(i).matrices(w,1) = array2table(sum(XOXO(i).matrices{w,:}));
    end 
end 

%STEP 7
clear Matrix_SumFlux_per_cluster_and_reactions_for_each_metabolites
for i= 1:numel(Sum_of_XOXO_columns) %width
    Matrix_SumFlux_per_cluster_and_reactions_for_each_metabolites(:,i) = Sum_of_XOXO_columns(i).matrices(:,1);
end 
Matrix_SumFlux_per_cluster_and_reactions_for_each_metabolites.Properties.RowNames(:)=PerfectMatrix.Properties.RowNames(:);
model.metNames(ismember(model.mets,'xolest2_hs[e]'))=cellstr('Cholesterol EsterFullR');
%changing the rownames to metabolite names with the help of simpleRev_recon3
pk = Matrix_SumFlux_per_cluster_and_reactions_for_each_metabolites.Properties.RowNames(:) ;
pk = strrep(pk,'[u]','[e]'); %--> replacing [u] by [e] due to the annotation in Recon3D being [e] instead of [u]
[I,IA,IB]=intersect(pk,model.mets);
Metabolites_names=pk;
Metabolites_names(IA) = model.metNames(IB);

Matrix_SumFlux_per_cluster_and_reactions_for_each_metabolites.Properties.RowNames(:)=Metabolites_names(:);
%changing the variable names to proper ones 
Matrix_SumFlux_per_cluster_and_reactions_for_each_metabolites.Properties.VariableNames(:) = Cluster_names_with_externalCompartment(:);

%STEP 8
%cleaning the matrix from 0 fluxes or irrelevant low fluxes otherwise the
%plotting is too complicated and you can't see anything

sum_of_matrix_rows = sum(table2array(Matrix_SumFlux_per_cluster_and_reactions_for_each_metabolites),2);
count = 0;
for i= 1:height(Matrix_SumFlux_per_cluster_and_reactions_for_each_metabolites)
    if sum_of_matrix_rows(i) > 1e-4
        count = count + 1;
        CleanedMatrix(count,:) = Matrix_SumFlux_per_cluster_and_reactions_for_each_metabolites(i,:);
        CleanedMatrix.Properties.RowNames(count) = Matrix_SumFlux_per_cluster_and_reactions_for_each_metabolites.Properties.RowNames(i);
    elseif sum_of_matrix_rows(i) < -1e-4
        count = count + 1;
        CleanedMatrix(count,:) = Matrix_SumFlux_per_cluster_and_reactions_for_each_metabolites(i,:);
        CleanedMatrix.Properties.RowNames(count) = Matrix_SumFlux_per_cluster_and_reactions_for_each_metabolites.Properties.RowNames(i);
    end 
end 
CleanedMatrix.Properties.VariableNames(:) = Cluster_names_with_externalCompartment(:);

%differential presence
Ranked_matrix = table(zeros(height(Matrix_SumFlux_per_cluster_and_reactions_for_each_metabolites),1));
for i = 1:height(Matrix_SumFlux_per_cluster_and_reactions_for_each_metabolites)
    maxValue = max(Matrix_SumFlux_per_cluster_and_reactions_for_each_metabolites{i,:});
    minValue = min(Matrix_SumFlux_per_cluster_and_reactions_for_each_metabolites{i,:});
if sortFlag==1
    Ranked_matrix(i,1) = table(maxValue - minValue);
elseif sortFlag==2  %TS alternative
        Ranked_matrix(i,1) = table(maxValue);
elseif sortFlag==3  %TS alternative
        Ranked_matrix(i,1) = table(minValue);
end
    Ranked_matrix.Properties.RowNames(i) = Matrix_SumFlux_per_cluster_and_reactions_for_each_metabolites.Properties.RowNames(i);
end 
if sortFlag==3
    Ranked_matrix_final = sortrows(Ranked_matrix,1,'ascend');
else
    Ranked_matrix_final = sortrows(Ranked_matrix,1,'descend');
end
inorganic_mets = ["Sodium","Water","Carbon Dioxide","Orthophosphate","Proton"];
Ranked_matrix_final(find(ismember(Ranked_matrix_final.Properties.RowNames(:),inorganic_mets)),:) = [];


%adding a column for the pathway of each reaction 
ModelIndxPathways = find(ismember(multicell_model.mets,PerfectMatrix.Properties.RowNames(:)));

clear Ranked_matrix_reactionANDpathways
Ranked_matrix_reactionANDpathways(:,1) = table(Ranked_matrix.Properties.RowNames);
Ranked_matrix_reactionANDpathways(:,2) = multicell_model.subSystems(ModelIndxPathways);
Ranked_matrix_reactionANDpathways(:,3) = table(Ranked_matrix.Var1);

Ranked_matrix_reactionANDpathways(find(ismember(Ranked_matrix_reactionANDpathways.Var1(:),inorganic_mets)),:) = [];
ColumnNames = ["Metabolites","Pathways","MaxFlux-MinFlux"];
Ranked_matrix_reactionANDpathways.Properties.VariableNames(:) = ColumnNames(:);
Ranked_matrix_reactionANDpathways = sortrows(Ranked_matrix_reactionANDpathways,3,'descend');

% writetable(Ranked_matrix_reactionANDpathways,'rankedmetabolitesWithPathways(FBA)','Delimiter','\t')
%taking the thrity first mets from the ranked final matrix
thirtyfivefirst_interesting_met = table(Ranked_matrix_final.Properties.RowNames(1:35,1)); %35
Interesting_met_matrix = Matrix_SumFlux_per_cluster_and_reactions_for_each_metabolites(find(ismember(table(Matrix_SumFlux_per_cluster_and_reactions_for_each_metabolites.Properties.RowNames),thirtyfivefirst_interesting_met)),:);

cgo_Path = clustergram(table2array(Interesting_met_matrix),...
        'RowLabels',Interesting_met_matrix.Properties.RowNames ,...
        'ColumnLabels', Interesting_met_matrix.Properties.VariableNames,...
        'ColumnLabelsRotate',340, ...
        'Cluster', 'all',...
        'LabelsWithMarkers',true,...
        'Colormap',redbluecmap,...
        'DisplayRange',10,...
        'DisplayRatio', [0.2 0.2]);

addTitle(cgo_Path,{'FBA of exchanged metabolites',['Model: ' num2str(modelNr) ', sortflag: ' num2str(sortFlag)]});
                                
         plot(cgo_Path);

         figureHandle = gcf;
         %# make all text in the figure to size 14 and bold
         fig_gcf = findall(0,'type','figure', 'tag', 'Clustergram');
         set(findall(figureHandle,'type','text'),'fontSize',16,'fontWeight','bold')
         set(findall(figureHandle, 'type','axes','tag','HeatMapAxes'),'fontsize',16)

% red (positive values) indicates the metabolite leaves the corresponding cluster
% blue (negative values) indicates the metabolite enters the corresponding cluster
data=Matrix_SumFlux_per_cluster_and_reactions_for_each_metabolites;
% data('Urea',:)
delete(['data' num2str(modelNr) '.xlsx'])
writetable(data,['data' num2str(modelNr) '.xlsx'],'WriteRowNames',1)
save(['FBAmodel' num2str(modelNr)],'FBAsolution','multicell_model')

%% check-out a specific metabolite
temp=find(ismember(multicell_model.mets,'glc_D[e]'));
temp=find(ismember(multicell_model.mets,'o2[e]'));
multicell_model.mets(temp)
multicell_model.S(temp,:);
temp2=find(multicell_model.S(temp,:));
printRxnFormula(multicell_model,multicell_model.rxns(temp2));
FBAsolution.v(temp2)

solf=FBAsolution.f;
resLactate=data;