clearvars -except solverOK
clc, close all force

O2uptake=21.253*5 %2; changed for REVISION
modelNrs=1:6
factors=[1]

resAll=[];
for counter2=1:numel(factors)
    factor=factors(counter2)
    res=[];
    for counter=1:numel(modelNrs)
        modelNr=modelNrs(counter)
        [solf,resEX] = FBA_mediumConc_varO2_v2_PAPER(modelNr,O2uptake*factor)
% % %         res(:,:,counter)=[repmat(solf,size(resEX,1),1), table2array(resEX)];
    end
    res
%     resAll(:,:,counter2)=res;
end
% resAll
try
    res=resEX({'Oxygen','Bicarbonate','D-Glucose','(S)-Lactate','Acetate','(R)-Mevalonate'},:)
end
try
    res=resEX({'Oxygen','Bicarbonate','D-Glucose','(S)-Lactate','(R)-3-Hydroxybutyrate','(R)-Mevalonate'},:)
end
try
    res=resEX({'Oxygen','Bicarbonate','D-Glucose','(S)-Lactate','(R)-3-Hydroxybutyrate','Acetate','(R)-Mevalonate'},:)
end
