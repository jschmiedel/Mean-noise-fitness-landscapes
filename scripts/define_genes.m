function define_genes()

%read gene IDs
GLU_reads = readtable('data/Keren2016_GSM2235555_CompetitionDataGlucose.txt','Delimiter','\t');
genes.ID = unique(cellfun(@(x) x(1:end-4),GLU_reads.Properties.VariableNames(3:end),'UniformOutput',false))';

%read gene wt promoter output from Keren2016
GLU_wtExpr = readtable('data/Keren2016_mmc3.xlsx','Sheet',2); %other sheets also have their binned data + fits
[a,b]=ismember(GLU_wtExpr.Gene,genes.ID);
genes.GLU_wtExpr(b(a),1) = GLU_wtExpr.wtExpressionGlucose(a);

save('data/genes.mat','genes')