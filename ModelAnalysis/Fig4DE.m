% Fig 4D: Predicted medium uptake and secretion rates of midbrain organoids for key metabolites.
clear all, close all, clc
[NUM,TXT,RAW]=xlsread('reEX_v2_factor1.xlsx');

%%
clc
O2=[NUM(1,2:7); NUM(11,2:7); NUM(22,2:7); NUM(32,2:7)]
Glc=[NUM(3,2:7); NUM(13,2:7); NUM(24,2:7); NUM(34,2:7)]
Lac=[NUM(4,2:7); NUM(14,2:7); NUM(25,2:7); NUM(35,2:7)]
HBut=[NUM(5,2:7); NUM(15,2:7); NUM(26,2:7); NUM(36,2:7)]
Ac=[NUM(6,2:7); NUM(16,2:7); NUM(27,2:7); NUM(37,2:7)]
Mev=[NUM(7,2:7); NUM(17,2:7); NUM(28,2:7); NUM(38,2:7)]

%% 4D
data=[Glc(:,6)'; O2(:,6)'; Lac(:,6)'; Ac(:,6)'; HBut(:,6)'; Mev(:,6)']
figure
bar(data)
legend('WT-D30','WT-D60','PD-D30','PD-D60','fontweight','bold','fontsize',10)
ylabel('medium secretion/uptake [a.u.]','fontweight','bold','fontsize',12)
ylim([-15 50])

set(gca,'xtick',1:6,'xticklabel',{'Glucose';'Oxygen';'Lactate';'Acetate';'3-Hydroxybutyrate';'(R)-Mevalonate'},'fontweight','bold','fontsize',12);
xtickangle(45)

%% 4E OLD
data=Lac'
figure
bar(data)
legend('WT-D30','WT-D60','PD-D30','PD-D60','fontweight','bold','fontsize',10)
% ylabel('lactate secretion/uptake [a.u.]','fontweight','bold','fontsize',12)
ylabel('lactate consumption/production [a.u.]','fontweight','bold','fontsize',12)

set(gca,'xtick',1:6,'xticklabel',{'GABAergic neurons';'Neurons';'Neural progenitors';'Dopaminergic neurons';'Astrocyte-like glia progenitors';'exchange with medium'},'fontweight','bold','fontsize',12);
xtickangle(45)

% + resize manually

%% 4E NEW
data=-Lac'
figure
bar(data)
legend('WT-D30','WT-D60','PD-D30','PD-D60','fontweight','bold','fontsize',10,'Location','southeast')
% ylabel('lactate secretion/uptake [a.u.]','fontweight','bold','fontsize',12)
ylabel('lactate production/consumption [a.u.]','fontweight','bold','fontsize',12)

set(gca,'xtick',1:6,'xticklabel',{'GABAergic neurons';'Neurons';'Neural progenitors';'Dopaminergic neurons';'Astrocyte-like glia progenitors';'exchange with medium'},'fontweight','bold','fontsize',12);
xtickangle(45)

% + resize manually

%% 4F
data=HBut'
figure
bar(data)
legend('WT-D30','WT-D60','PD-D30','PD-D60','fontweight','bold','fontsize',10)
ylabel({'3-hydroxybutyrate';'consumption/production [a.u.]'},'fontweight','bold','fontsize',12)

set(gca,'xtick',1:6,'xticklabel',{'GABAergic neurons';'Neurons';'Neural progenitors';'Dopaminergic neurons';'Astrocyte-like glia progenitors';'exchange with medium'},'fontweight','bold','fontsize',12);
xtickangle(45)

% + resize manually