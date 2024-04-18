% driver analysis, TS 07/23
%% scFASTCORMICS HPC results Maria Pacheco, 10.07.23
clearvars -except solverOK
clc, close all
load('Claudia_medium_cons_250723.mat')
Results_keep

%% model sizes
res=[];
for counter=1:numel(Results_keep)
%     model=Results_keep(counter).multi_cell_population_model;
    model=Results_keep(counter).multi_cell_population;
    res=[res; counter numel(model.rxns)];
end
disp('Model sizes:')
disp(res)

%% model sizes per cluster (reactions)
res=[];
for counter=1:numel(Results_keep)
%     model=Results_keep(counter).multi_cell_population_model;
    model=Results_keep(counter).multi_cell_population;
    temp=[];
    for counter2=1:5
        idx = strfind(model.rxns, ['_' num2str(counter2)]);
        idx = find(not(cellfun('isempty', idx)));
        temp=[temp, numel(idx)];
    end
    res=[res; temp];
end
disp('Model sizes per clusters (cols):')
disp(res)
namesCluster={'c1','c2','c3','c45','c6'}

min(min(res(1:4,:)))
max(max(res(1:4,:)))

%% model sizes per cluster (metabolites)
res=[];
for counter=1:numel(Results_keep)
%     model=Results_keep(counter).multi_cell_population_model;
    model=Results_keep(counter).multi_cell_population;
    temp=[];
    for counter2=1:5
        idx = strfind(model.mets, ['_' num2str(counter2)]);
        idx = find(not(cellfun('isempty', idx)));
        temp=[temp, numel(idx)];
    end
    res=[res; temp];
end
disp('Model sizes per clusters (cols):')
disp(res)
namesCluster={'c1','c2','c3','c45','c6'}

min(min(res(1:4,:)))
max(max(res(1:4,:)))

%% Jaccard plot Whole Models (FOR PAPER, only 4 conditions)
res=nan(4,4); %nan(numel(Results_keep));
for counter=1:4 %numel(Results_keep)
    for counter2=1:4 %numel(Results_keep)
        A1=find(Results_keep(counter).A);
        A2=find(Results_keep(counter2).A);
        res(counter,counter2)=numel(intersect(A1,A2))/numel(union(A1,A2));
    end
end
% disp('Jaccard similarity Whole Models:')
disp(res)

J=res;
names_col={'WT-D30';'WT-D60'; 'PD-D30'; 'PD-D60'} %{'Ctrl30';'Ctrl60'; 'PD30'; 'PD60'};
altcolor = [255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;...
    255 0 0; 204 0 0; 152 0 0; 102 0 0; 51 0 0]/255;
cgo_J = clustergram(J,...
    'RowLabels', names_col,...
    'ColumnLabels', names_col,...
    'ColumnLabelsRotate',45, ...
    'Cluster', 'all', ...
    'Annotate', 'true',...
    'symmetric','False',...
    'AnnotColor','k',...
    'Colormap', altcolor)

% addTitle(cgo_J,{'Jaccard similarity Whole Models'});
plot(cgo_J);

% clear gcf
figureHandle = gcf;
%# make all text in the figure to size 14 and bold
set(findall(figureHandle,'type','text'),'fontSize',18,'fontWeight','bold')
% gcf.Position=[400 300 700 500]

%and manually change label fontsize to 12 and font to bold
