function promoter_expression_Sharon2014()

load data/promoters.mat
load data/genes.mat

Sharon2014_promoters = readtable('data/Sharon2014_TableS7.txt','Delimiter','\t');
%mean and noise (CV^2 !!!) when grown on glucose on linear scale

%map IDs
[a,b]=ismember(Sharon2014_promoters.LibID,promoters.ID);

%average over replicates
avg_mean = log2(mean([Sharon2014_promoters.Replicate1_ExpressionMean Sharon2014_promoters.Replicate2_ExpressionMean],2));
avg_noise = log2(sqrt(mean([Sharon2014_promoters.Replicate1_ExpressionNoise Sharon2014_promoters.Replicate2_ExpressionNoise],2)));

%error between replicates
error_mean = noise([Sharon2014_promoters.Replicate1_ExpressionMean Sharon2014_promoters.Replicate2_ExpressionMean],2,'cv');
error_noise = noise(sqrt([Sharon2014_promoters.Replicate1_ExpressionNoise Sharon2014_promoters.Replicate2_ExpressionNoise]),2,'cv');
%calculate expected error from read counts
error_reads = sqrt((1./Sharon2014_promoters.Replicate1_NumReads) + (1./Sharon2014_promoters.Replicate2_NumReads)); %checked

%model (avg) error as function of counting error
error_mean_modeled = smooth(error_reads,error_mean,0.5,'loess');
error_noise_modeled = smooth(error_reads,error_noise,0.5,'loess');

%adjust mean output between Sharon2014 and Keren2016
P = polyfit(avg_mean(a),promoters.GLU_expression_Keren2016,1);
avg_mean_adj = avg_mean*P(1) + P(2);

%write data
promoter_expression = table();
promoter_expression.mean(b(a),1) = avg_mean_adj(a);
promoter_expression.noise(b(a),1) = avg_noise(a);
promoter_expression.error_mean(b(a),1) = error_mean(a);
promoter_expression.error_noise(b(a),1) = error_noise(a);
promoter_expression.error_mean_modeled(b(a),1) = error_mean_modeled(a);
promoter_expression.error_noise_modeled(b(a),1) = error_noise_modeled(a);

promoter_expression.Properties.Description = 'promoter_expression_GLU_Sharon2014';

save('data/promoter_expression_GLU_Sharon2014.mat','promoter_expression')


%promoters that differ too much between Keren2016 and Sharon2014
g0p5 = (abs((avg_mean_adj(a)-promoters.GLU_expression_Keren2016))>0.5);

f=figure;
%plot mean & noise of promoters used in Keren2016 versus all used in Sharon2014
subplot(2,2,1)
plot(avg_mean,avg_noise,'.')
hold on
plot(avg_mean(a),avg_noise(a),'r.','MarkerSize',15)
xlabel('mean output (log2)')
ylabel('noise (CV) (log2)')
title('promoter mean:noise relationship')
legend('all promoters','Keren2016','Location','SouthWest')
box off

%compare Sharon2014 to Keren2016 scale
s=subplot(2,2,2);
avg_S = avg_mean(a);
plot(promoters.GLU_expression_Keren2016,avg_S,'k.','MarkerSize',12)
lsline
hold on
plot(promoters.GLU_expression_Keren2016(g0p5),avg_S(g0p5),'r.','MarkerSize',12)
title('expression strength comparision (log2)')
xlabel('Keren2016')
ylabel('Sharon2014')
[r,p] = corr(avg_S,promoters.GLU_expression_Keren2016,'type','Pearson');
text(0,4,['all, R^2 = ' num2str(r^2,2) ', p = ' num2str(p,1)],'FontSize',9)
[r,p] = corr(avg_S(~g0p5),promoters.GLU_expression_Keren2016(~g0p5),'type','Pearson');
text(0,3.3,['w/o outliers, R^2 = ' num2str(r^2,2) ', p = ' num2str(p,1)],'FontSize',9)
set(s,'XTick',0:2:6,'YTick',-2:2:4,'FontSize',10)
box off
axis([-0.5 7 -3.5 4.5])

%analyse error in mean
subplot(2,2,3)
loglog(error_reads,error_mean,'.')
hold on
loglog(error_reads(a),error_mean(a),'r.')
plot(error_reads,error_mean_modeled,'k.')
title('mean expression')
xlabel('counting error|read counts')
ylabel('CV of mean')
axis([0 inf 0 inf])
legend('all prom.','Keren2016 prom.','modeled error','Location','NorthWest')
box off

%analyse error in noise
subplot(2,2,4)
loglog(error_reads,error_noise,'.')
hold on
loglog(error_reads(a),error_noise(a),'r.')
plot(error_reads,error_noise_modeled,'k.')
title('noise')
xlabel('counting error|read counts')
ylabel('CV of noise')
axis([0 inf 0 inf])
box off
%
print(f,'-depsc2','results/setup_Sharon2014_promoter_expression_GLU.eps')

