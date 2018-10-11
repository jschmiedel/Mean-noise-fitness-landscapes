%% assemble table with experimental data on expression-senstivity and noise of genes
function assemble_largescale_data_for_validation()
load data/genes.mat

%ORF - gene name conversion
ORF2gene = readtable('data/yeast_gene_ORF2.csv','TreatAsEmpty',{'NA','',' '},'ReadVariableNames',false);
ORF2gene = ORF2gene(~strcmp('',ORF2gene.Var2),1:2);
ORF2gene.Properties.VariableNames = {'ORF','gene'};

for i=1:numel(genes.ID)
    match = find(strcmp(genes.ID{i},ORF2gene.gene));
    if ~isempty(match)
        conv_Keren_SGD_idx(i,1) = match;
    else
        if strcmp(genes.ID{i},'TIM23')
            conv_Keren_SGD_idx(i,1) = find(strcmp('MAS6',ORF2gene.gene));
        elseif strcmp(genes.ID{i},'TIM10')
            conv_Keren_SGD_idx(i,1) = find(strcmp('MRS11',ORF2gene.gene));
        elseif strcmp(genes.ID{i},'USV1')
            conv_Keren_SGD_idx(i,1) = find(strcmp('YPL230W',ORF2gene.gene));
        end  
    end
end
endogenous_validation = ORF2gene(conv_Keren_SGD_idx,:);

%essential genes
output = readtable('data/essential_ORFs.txt','Delimiter','\t','TreatAsEmpty',{'NA','',' '});
essential_ORFs = strtrim(output.ORF_name(2:end));
[a,b] = ismember(essential_ORFs,endogenous_validation.ORF);
endogenous_validation = [endogenous_validation array2table(false(size(endogenous_validation,1),1),'VariableNames',{'essential'})];
endogenous_validation.essential(b(a)) = true;

%over-expression sensitive genes
endogenous_validation = [endogenous_validation array2table(false(size(endogenous_validation,1),1),'VariableNames',{'OE_sensitive'})];
load('data/Sopko2006_OEsensitive.mat')
[a,b] = ismember(Sopko2006.ORF,endogenous_validation.ORF);
endogenous_validation.OE_sensitive(b(a)) = true;

load('data/Makanae2013_OEsens.mat')
[a,b] = ismember(Makanae2013.ORF,endogenous_validation.ORF);
endogenous_validation.OE_sensitive(b(a)) = true;

%noise data from Newman2006
output = readtable('data/Newman2006_TableS1_nature04785-s03.txt','Delimiter','\t','Headerlines',2,'TreatAsEmpty',{'NA','',' '});
endogenous_validation = [endogenous_validation array2table(NaN(size(endogenous_validation,1),1),'VariableNames',{'Noise_Newman_DM_SD'})];
endogenous_validation = [endogenous_validation array2table(NaN(size(endogenous_validation,1),1),'VariableNames',{'Noise_Newman_DM_YPD'})];
[a,b] = ismember(output.ORF,endogenous_validation.ORF);
endogenous_validation.Noise_Newman_DM_SD(b(a)) = output.DMSD(a);
endogenous_validation.Noise_Newman_DM_YPD(b(a)) = output.DMYEPD(a);

%noise data from Stewart-Ornstein 2012
endogenous_validation = [endogenous_validation array2table(NaN(size(endogenous_validation,1),1),'VariableNames',{'Noise_StewartO_DevTot'})];
load('data/StewardOrnstein2012_Noise.mat')
% calculate DM measure
Nisnan = ~isnan(StewartO2012.ExprNoise.P_mean) & ~isnan(StewartO2012.ExprNoise.CvTot);
smooth_span = 30;
Q = .5;
[smooth0,~,sort_i] = smooth_quantiles(StewartO2012.ExprNoise.P_mean(Nisnan),StewartO2012.ExprNoise.CvTot(Nisnan),smooth_span,Q,1);
clear a
a(sort_i) = smooth0;
StewartO2012.ExprNoise.DevTot = StewartO2012.ExprNoise.CvTot(Nisnan)./a';
StewartO2012.ExprNoise.DevTot(StewartO2012.ExprNoise.DevTot==0) = NaN;
StewartO2012.ExprNoise.DevTot(StewartO2012.ExprNoise.DevTot==inf) = NaN;
[a,b] = ismember(StewartO2012.ORF,endogenous_validation.ORF);
endogenous_validation.Noise_StewartO_DevTot(b(a)) = StewartO2012.ExprNoise.DevTot(a);

endogenous_validation.Properties.RowNames = genes.ID;

save('data/endogenous_validation_data.mat','endogenous_validation');