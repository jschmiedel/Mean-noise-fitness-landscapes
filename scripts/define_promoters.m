function define_promoters()

%promoter IDs
GLU_reads = readtable('data/Keren2016_GSM2235555_CompetitionDataGlucose.txt','Delimiter','\t');
[promoters.ID,idx] = sort(GLU_reads.PromoterID_SharonEtAl__2012_,'ascend');

%promoter mean output as estimated in Keren2016
%glucose
promoters.GLU_expression_Keren2016 = GLU_reads.PromoterExpression_log2_(idx);

%promoter description and sequence from Sharon2012
Sharon2012_PromSeqs0 = readtable('data/Sharon2012_TableS3.xlsx'); %% copy-pasted from Sharon2012 supplement
[a,b]=ismember(Sharon2012_PromSeqs0.LibraryID,promoters.ID);
promoters.Description(b(a),1) = Sharon2012_PromSeqs0.Description(a);
promoters.OligoSequence(b(a),1) = Sharon2012_PromSeqs0.OligoSequence(a);

save('data/promoters.mat','promoters')