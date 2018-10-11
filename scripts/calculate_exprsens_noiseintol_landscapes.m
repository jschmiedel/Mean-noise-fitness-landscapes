function calculate_exprsens_noiseintol_landscapes()
% calculate expression-sensitivity and noise-intolerance from smoothed landscapes
% as derivates of expression and noise around wild-type

load data/genes.mat
load data/promoters.mat
load data/exclusion_gene_promoters.mat

load data/promoter_expression_GLU_Sharon2014.mat
load data/fitness_t123_glucose.mat

load data/FNM_smooth_fun_grid05_025.mat
load data/FNM_smooth_fun_grid05_025_landscape_PCA.mat

Fitness = fitness{:,:};

show_genes = genes.ID(genes.GLU_wtExpr > 3 & genes.GLU_wtExpr < 5);

%randomizations of landscapes
nr_randomizations = 10000;
%randomizatoin vectors
rng(1603)
zr = arrayfun(@(x) randperm(sum(sum(full(exclusion_gene_promoters(:,ismember(genes.ID,show_genes)))==0,2)>0)),1:nr_randomizations,'uni',false);

%% noise intolerance: sensitivity to noise at wild-type expression
MODE.local_grid_n = 21;

%first for noise intolerance
noise_intolerance_real = NaN(1,numel(show_genes));
noise_intolerance_rand = NaN(nr_randomizations,numel(show_genes));

noise_intolerance_real2 = NaN(1,numel(show_genes));
noise_intolerance_rand2 = NaN(nr_randomizations,numel(show_genes));

noise_intolerance_real3 = NaN(1,numel(show_genes));
noise_intolerance_rand3 = NaN(nr_randomizations,numel(show_genes));

for i=1:numel(show_genes)
    
    gene_idx = find(strcmp(genes.ID,show_genes(i)));
    MODE.xv = genes.GLU_wtExpr(gene_idx);
    MODE.n = numel(MODE.xv);
    MODE.yv = -3 : .05 : -1;
    
    % smooth landscapes with gaussian blurring and numerically estimate gradients
    %exclude data points
    Nisnan = exclusion_gene_promoters(:,gene_idx)==0 & ~isnan(promoter_expression.mean);
    Nisnan = Nisnan & promoter_expression.mean>2 & promoter_expression.mean<6;
    
    %define data points
    % as averages
    x=promoter_expression.mean(Nisnan);
    y=promoter_expression.noise(Nisnan);
    z=Fitness(Nisnan,gene_idx);
    err_x = promoter_expression.error_mean_modeled(Nisnan);
    err_y = promoter_expression.error_noise_modeled(Nisnan);
    err_z = error_fitness{Nisnan,gene_idx};
    
    %calculate fitness on stripe along noise at wild-type expression
    F = fitness_noise_mean_smoothing_v5(x,y,z,err_x,err_y,err_z,MODE);
    
    
    %noise intolerance as maximal -slope
    noise_intolerance_real(1,i) = nanmax(-gradient(F,diff(MODE.yv([1 2]))));
    
    %calculate noise intolerance as mean of -(piece-wise slopes)
    noise_intolerance_real2(1,i) = nanmean(-gradient(F,diff(MODE.yv([1 2]))));
    
    %mean noise curvature
    noise_intolerance_real3(1,i) = nanmean(-gradient(gradient(F,diff(MODE.yv([1 2]))),diff(MODE.yv([1 2]))));
    
    
    %now randomize fitness values to get background distribution
    for j=1:nr_randomizations
        if mod(j,1000) == 0
            sprintf('%0.1f%%',100*(((i-1)*nr_randomizations + j)/(numel(show_genes)*nr_randomizations)))
        end
        
        z_rand = z(zr{j}(zr{j} <= sum(Nisnan)));
        err_z_rand = err_z(zr{j}(zr{j} <= sum(Nisnan)));
        
        %calculate fitness on stripe along noise at wild-type expression
        F = fitness_noise_mean_smoothing_v5(x,y,z_rand,err_x,err_y,err_z_rand,MODE);
    
        %calculate noise intolerance as mean of -(piece-wise slopes)
        noise_intolerance_rand(j,i) = nanmax(-gradient(F,diff(MODE.yv([1 2]))));
        noise_intolerance_rand2(j,i) = nanmean(-gradient(F,diff(MODE.yv([1 2]))));
        noise_intolerance_rand3(j,i) = nanmean(-gradient(gradient(F,diff(MODE.yv([1 2]))),diff(MODE.yv([1 2]))));
    end
end

%% expression-sensitivity: sensitivity to mean expression at diff. noise levels
expression_sensitivity_real = NaN(1,numel(show_genes));
expression_sensitivity_rand = NaN(nr_randomizations,numel(show_genes));

for i=1:numel(show_genes)
    
    gene_idx = find(strcmp(genes.ID,show_genes(i)));
    %for each gene pre-calculate weights on grid/stripe
    MODE.xv = (genes.GLU_wtExpr(gene_idx)-1.5) : 0.05 : (genes.GLU_wtExpr(gene_idx)+1.5);
    MODE.n = numel(MODE.xv);
    MODE.yv = -3;
    
    % smooth landscapes with gaussian blurring and numerically estimate gradients
    %exclude data points
    Nisnan = exclusion_gene_promoters(:,gene_idx)==0 & ~isnan(promoter_expression.mean);
    Nisnan = Nisnan & promoter_expression.mean>2 & promoter_expression.mean<6;
    
    %define data points
    % as averages
    x=promoter_expression.mean(Nisnan);
    y=promoter_expression.noise(Nisnan);
    z=Fitness(Nisnan,gene_idx);
    err_x = promoter_expression.error_mean_modeled(Nisnan);
    err_y = promoter_expression.error_noise_modeled(Nisnan);
    err_z = error_fitness{Nisnan,gene_idx};
    
    %calculate fitness on stripe along noise at wild-type expression
    F = fitness_noise_mean_smoothing_v5(x,y,z,err_x,err_y,err_z,MODE);
    ES_d2curve(:,i) = -gradient(gradient(F(1,:),diff(MODE.xv([1 2]))),diff(MODE.xv([1 2])));
    
    % expression sensitivity as mean curvature
    expression_sensitivity_real(1,i) = nanmean(ES_d2curve(:,i));

    %now randomize fitness values to get background distribution
    for j=1:nr_randomizations
        if mod(j,1000) == 0
            sprintf('%0.1f%%',100*(((i-1)*nr_randomizations + j)/(numel(show_genes)*nr_randomizations)))
        end
        
        z_rand = z(zr{j}(zr{j} <= sum(Nisnan)));
        err_z_rand = err_z(zr{j}(zr{j} <= sum(Nisnan)));
        
        %calculate fitness on stripe along noise at wild-type expression
        F0 = fitness_noise_mean_smoothing_v5(x,y,z_rand,err_x,err_y,err_z_rand,MODE);
    
        helper = -gradient(gradient(F0(1,:),diff(MODE.xv([1 2]))),diff(MODE.xv([1 2])));
        expression_sensitivity_rand(j,i) = nanmean(helper);
    end
end

save('data/FNM_smooth_fun_grid05_025_slice_10000randomizations_sensitivity_analysis.mat',...
    'noise_intolerance_real','noise_intolerance_real2','noise_intolerance_real3',...
    'noise_intolerance_rand','noise_intolerance_rand2','noise_intolerance_rand3',...
    'expression_sensitivity_real','expression_sensitivity_rand')

%% load data
load data/FNM_smooth_fun_grid05_025_slice_10000randomizations_sensitivity_analysis.mat

%% example calculation of expression-sensitivity and noise-intolerance
F = principal_landscapes(:,:,1) + principal_landscapes(:,:,2);

Fxx = -gradient(gradient(F(1,:),0.05),0.05);
Fy = -gradient(F(:,31),0.025);
xv = -1.5:.05:1.5;
yv = -3:.025:-1;

xb = [-1 1];

rc=4;

f=figure;
s=subplot(rc,rc,1);
plot(xv,F(1,:))
xlabel('mean expression')
ylabel('fitness')
box off
set(s,'XTick',-1:1:1,'YTick',-0.04:0.02:0)
axis([-1.5 1.5 -0.04 0.005])

s=subplot(rc,rc,rc+1);
plot(xv(3:end-2),Fxx(3:end-2))
hold on
plot([xv(3) xv(end-2)],ones(2,1)*mean(Fxx(3:end-2)),'r')
xlabel('mean expression')
ylabel('\partial^2fitness/\partialmean^2')
box off
set(s,'XTick',[-1:1:1],'YTick',0:0.01:0.03)
axis([-1.5 1.5 -0.005 0.025])

s=subplot(rc,rc,2);
plot(yv',F(:,31))
xlabel('noise')
ylabel('fitness')
box off
set(s,'XTick',[-3:1:-1],'YTick',-0.04:0.02:0)
axis([-3 -1 -0.04 0.005])

s=subplot(rc,rc,rc+2);
plot(yv,Fy)
hold on
plot(yv,ones(numel(yv),1)*max(Fy),'r')
xlabel('noise')
ylabel('\partialfitness/\partialnoise')
set(s,'XTick',[-3:1:-1],'YTick',0:0.01:0.03)
box off
axis([-3 -1 -0.005 0.035])

print(f,'-depsc2','results/ESandNI_sketch.eps')

%% compare noise_intolerance metrics

f=figure;
s=subplot(2,2,1);
scatter(noise_intolerance_real,noise_intolerance_real2,'filled')
axis tight
xlabel('maximal negative slope')
ylabel('average negative slope')
set(s,'XTick',0:0.01:0.05,'YTick',-0.01:0.01:0.05,'FontSize',10)

s=subplot(2,2,2);
scatter(noise_intolerance_real,noise_intolerance_real3,'filled')
axis tight
xlabel('maximal negative slope')
ylabel('average negative curvature')
set(s,'XTick',0:0.01:0.05,'YTick',-0.01:0.005:0.05,'FontSize',10)

print(f,'-depsc2','results/comparison_NImetrics.eps')

%% compare noise intolerance and expression-sensitivity between genes

NI = noise_intolerance_real; %maximal slope
NI_rand = noise_intolerance_rand;
%calculate pval and fdr
NI_pval_genewise = arrayfun(@(x) sum(NI_rand(:,x) >= NI(1,x))/(nr_randomizations)...
    ,1:numel(show_genes));
NI_fdr_genewise = mafdr(NI_pval_genewise,'BHFDR',true);

ES = expression_sensitivity_real; %mean curvature without boundary regions
ES_rand = expression_sensitivity_rand;
%calculate pval and fdr
ES_pval_genewise = arrayfun(@(x) sum(ES_rand(:,x) >= ES(1,x))/(nr_randomizations)...
    ,1:numel(show_genes));
ES_fdr_genewise = mafdr(ES_pval_genewise,'BHFDR',true);

%plot xylim
xylim = [-.0025 0.031 -.01 .055];

f=figure;
subplot(4,4,[5 15])
h=gscatter(ES,NI,1:33,'k','o',10,'off','filled');
hold on
for i=1:numel(h)
    %mark FDR val on genes
    if NI_fdr_genewise(i) < 0.1 && ES_fdr_genewise(i) < 0.1
        marker = 'd';
    elseif NI_fdr_genewise(i) < 0.1 && ES_fdr_genewise(i) > 0.1
        marker = '^';
    elseif NI_fdr_genewise(i) > 0.1 && ES_fdr_genewise(i) < 0.1
        marker = '>';
    else 
        marker = 'o';
    end
    set(h(i),'MarkerFaceColor',[.5 .5 .5],'Marker',marker);
end

%plot coordinate system
plot(xylim([1 2]),zeros(2,1),'k--')
plot(zeros(2,1),xylim([3 4]),'k--')

%correlation between ES and NI
r = corr(ES', NI','type','Pearson');
%correlation in randomized data
r_rand = cell2mat(arrayfun(@(x) corr(ES_rand(x,:)',NI_rand(x,:)','type','Pearson'),1:nr_randomizations,'uni',false));
text(.02,.03,['R = ' num2str(r,2), ', p = ' num2str(1-sum(r_rand<r)/nr_randomizations,2)])

axis(xylim)
xlabel('expression sensitivity','FontSize',12)
ylabel('noise intolerance','FontSize',12)
set(gca,'XTick',0:0.01:xylim(2),'YTick',0:0.01:xylim(4),'FontSize',12)
box off

%expression-sensitivity histogram
s=subplot(4,4,[1 3]);
[hfdr,~]=hist(ES(ES_fdr_genewise < 0.1),-0.00625:0.0025:0.045);
[hnfdr,x]=hist(ES(ES_fdr_genewise > 0.1),-0.00625:0.0025:0.045);
bar(x,[hfdr; hnfdr]',1,'stacked')
hold on
set(s,'XTick',[],'YTick',0:5:15,'FontSize',12)
axis([xylim([1 2]) 0 12])
box off
text(0.005,11,[num2str(sum(ES_fdr_genewise < .1)) '/' num2str(numel(show_genes)) ' genes'])

%noise-tolerance histogram
s=subplot(4,4,[8 16]);
[hfdr,~]=hist(NI(NI_fdr_genewise < 0.1),-0.0325:0.005:0.075);
[hnfdr,x]=hist(NI(NI_fdr_genewise > 0.1),-0.0325:0.005:0.075);
barh(x,[hfdr; hnfdr]',1,'stacked')
hold on
set(s,'YTick',[],'XTick',0:5:15,'FontSize',12)
axis([ 0 13  xylim([3 4])])
box off
text(10,0.01,[num2str(sum(NI_fdr_genewise < .1)) '/' num2str(numel(show_genes)) ' genes'])

print(f,'-depsc2','results/ESvsNI.eps')


%% noise intolerance and expressionsensitivity from raw data & partial correlations

for i=1:numel(show_genes)
    gene = find(strcmp(genes.ID,show_genes{i}));
    Nisnan = exclusion_gene_promoters(:,gene)==0& promoter_expression.mean > genes.GLU_wtExpr(gene)-1.5 & promoter_expression.mean < genes.GLU_wtExpr(gene)+1.5;
    [r,p]=partialcorr([promoter_expression.mean(Nisnan) promoter_expression.noise(Nisnan) Fitness(Nisnan,gene)],'type','Pearson');
    FNM_pcorr(i,[1 2]) = r([1 2],3);
    FNM_pcorr_pval(i,[1 2]) = p([1 2],3);
end

FNM_pcorr(:,2) =  -FNM_pcorr(:,2);
%fdr from pvals
FDR_pc(:,1) = mafdr(FNM_pcorr_pval(:,1),'BHFDR',true);
FDR_pc(:,2) = mafdr(FNM_pcorr_pval(:,2),'BHFDR',true);
%mark genes significant at 10% FDR
sig = repmat({'Nsig'},size(FDR_pc));
sig(FDR_pc(:,1) < 0.1,1) = repmat({'sig'},sum(FDR_pc(:,1) < 0.1),1);
sig(FDR_pc(:,2) < 0.1,2) = repmat({'sig'},sum(FDR_pc(:,2) < 0.1),1);

f=figure;
subplot(4,4,[5 15])
h=scatter(FNM_pcorr(:,1),FNM_pcorr(:,2),50,'filled');
for i=1:numel(h)
    set(h(i),'MarkerEdgeColor','k','MarkerFaceColor','k','Marker','o');
end
hold on
plot([-1 1],zeros(2,1),'k--')
plot(zeros(2,1),[-1 1],'k--')

plot([-1 1],ones(2,1)*0.25,':','Color',ones(3,1)*0.6)
text(.6,.27,'FDR < 0.1','FontSize',12)

plot(ones(2,1)*0.25,[-1 1],':','Color',ones(3,1)*0.6)
text(.26,.7,'FDR < 0.1','FontSize',12)

[r,p]=corr(FNM_pcorr(:,1),FNM_pcorr(:,2),'type','Pearson');
text(.5,.5,['R = ' num2str(r,2), ', p = ' num2str(p,1)])

xlabel('expression sensitivity \rho_{(fitness,mean) \cdot noise}','FontSize',20)
ylabel('noise intolerance -\rho_{(fitness,noise) \cdot mean}','FontSize',20)
set(gca,'XTick',0:.25:.75,'YTick',0:.25:0.75,'FontSize',16)
axis([-.2 .75 -.2 .75])
box off

%expression-sensitivity histogram
s=subplot(4,4,[1 3]);
[h,x]=hist(FNM_pcorr(:,1),[-.0625:.125:.625]);
bar(x,h,1,'FaceColor',ones(3,1)*0.5)
hold on
plot(ones(2,1)*0.25,[0 12],':','Color',ones(3,1)*0.6)
set(s,'YTick',[],'XTick',0:5:15,'FontSize',12)
axis([-0.2 .75 0 12  ])
box off
text(.5,10,[num2str(sum(FNM_pcorr(:,1) > .25)) '/' num2str(size(FNM_pcorr,1)) ' genes'])

%noise-tolerance histogram
s=subplot(4,4,[8 16]);
[h,x]=hist(FNM_pcorr(:,2),[-.0625:.125:.5625]);
barh(x,h,1,'FaceColor',ones(3,1)*0.5)
hold on
plot([0 12],ones(2,1)*0.25,':','Color',ones(3,1)*0.6)
set(s,'XTick',[],'YTick',0:5:15,'FontSize',12)
axis([ 0 12  -.2 0.75])
box off
text(3,.65,[num2str(sum(FNM_pcorr(:,2) > .25)) '/' num2str(size(FNM_pcorr,1)) ' genes'])

print(f,'-depsc2',['results/raw_data_partialcorr_colored.eps'])

%% how well do principal topology loadings explain noise intolerance?
lm_results = LinearModel.fit(principal_landscapes_genescores(:,[1 2]),NI');
[R,p] = corr([NI',...
    lm_results.Coefficients.Estimate(2)*principal_landscapes_genescores(:,1)+ ...
    lm_results.Coefficients.Estimate(3)*principal_landscapes_genescores(:,2)],...
    'type','Pearson')
R.^2

f=figure;
h=lm_results.plot;
set(h(1),'MarkerFaceColor','k','MarkerEdgeColor','k','Marker','o')
axis tight
box off
ylabel('noise intolerance')
xlabel(['y ~ ',num2str(lm_results.Coefficients.Estimate(1),2), ' + ',num2str(lm_results.Coefficients.Estimate(2),2),'*PT1 + ',num2str(lm_results.Coefficients.Estimate(3),2),'*PT2' ])
title('')
set(gca,'xtick',[0:0.01:0.05]/0.0185842,'xticklabel',0:0.01:0.05)

print(f,'-depsc2','results/NI_PTreconstructed.eps')

%% compare measures against classifications from large scale datasets

load data/endogenous_validation_data.mat
idx_i = ismember(endogenous_validation.Properties.RowNames,genes.ID(genes.GLU_wtExpr > 3 & genes.GLU_wtExpr < 5));
endogenous_validation_limited = endogenous_validation(idx_i,:);

%%% and Keren2016 dosage-sensitivity
%read wild-type expression
GLU_wtExpr = readtable('data/Keren2016_mmc3.xlsx','Sheet',2); %other sheets also have their binned data + fits
GLU_fitness_fit_Keren2016 = readtable('data/Keren2016_mmc3.xlsx','Sheet',5,'ReadVariableNames',false);
Keren2016_binned_expr = GLU_fitness_fit_Keren2016{1,2:end};
helper_genes = GLU_fitness_fit_Keren2016{2:end,1};
Keren2016_binned_fitnessFIT = GLU_fitness_fit_Keren2016{2:end,2:end};
% in main paper fitness threshold is 5%
% fitness_thresholds = [0.02 0.05 0.1]; %2% and 10% give same results >> only show 5%
fitness_thresholds = 0.05;
clear Keren2016_wtFitness Keren2016_expr_close2wt Keren2016_dosageSensitivity Keren2016_sensitive_expr
for i=1:numel(show_genes)

    I = find(strcmp(helper_genes,show_genes(i)));
    %for each gene, find expression bin closest to wild-type expression
    [~,Keren2016_expr_close2wt] = min(abs(Keren2016_binned_expr - GLU_wtExpr.wtExpressionGlucose(I)));

    %find bins with less than 5% of wild-type fitness
    fitness_below_5p = find((Keren2016_binned_fitnessFIT(I,Keren2016_expr_close2wt) - Keren2016_binned_fitnessFIT(I,:)) > fitness_thresholds);

    %find bin w/ less than 5% fitness closests to wild-type
    if ~isempty(fitness_below_5p)
        [~,idx] = min(abs(fitness_below_5p - Keren2016_expr_close2wt));
        Keren2016_dosageSensitivity(i,1) = abs(Keren2016_binned_expr(idx) - GLU_wtExpr.wtExpressionGlucose(I));
        Keren2016_sensitive_expr(i,1) = Keren2016_binned_expr(idx);
    else
        Keren2016_dosageSensitivity(i,1) = inf;
        Keren2016_sensitive_expr(i,1) = NaN;
    end
end

%% compare expression-sensitivity, noise-tolerance and principal topology loadings with known dosage-sensitivities

matrix_es = repmat(ES',1,2);
matrix_nt = repmat(NI',1,2);

matrix_pt1 = repmat(principal_landscapes_genescores(:,1),1,2);
matrix_pt2 = repmat(principal_landscapes_genescores(:,2),1,2);

ess = sum(endogenous_validation_limited{:,[3]},2)>0;
oe = sum(endogenous_validation_limited{:,[4]},2)>0;
ess_oe = sum(endogenous_validation_limited{:,[3 4]},2)>0;
matrix_es(ess_oe,1) = NaN;
matrix_nt(ess_oe,1) = NaN;
matrix_es(~ess_oe,2) = NaN;
matrix_nt(~ess_oe,2) = NaN;

matrix_pt1(ess,1) = NaN;
matrix_pt2(oe,1) = NaN;
matrix_pt1(~ess,2) = NaN;
matrix_pt2(~oe,2) = NaN;

f=figure;
s=subplot(1,2,1);
h=plotSpread(matrix_es);
set(h{1},'MarkerSize',15,'Color','k')
hold on
boxplot(matrix_es)
h=findobj(gca,'tag','Outliers'); % Get handles for outlier lines.
set(h,'Marker','none')
h=refline(0,0);
set(h,'LineStyle','--','Color','k')

%calculate AUCs and pvals
[~,~,~,AUC] = perfcurve(ess_oe,ES,true);
p = ranksum(matrix_es(:,1),matrix_es(:,2),'tail','left');
text(.8,0.028,['AUC = ' num2str(AUC,2) ', p = ' num2str(p,2)],'FontSize',14)

set(s,'XTick',[1 2],'XTickLabel',{['insensitive, n=' num2str(sum(~ess_oe))],['essential | OE sens., n=' num2str(sum(ess_oe))]},'XTickLabelRotation',45,...
    'YTick',0:.01:.05,'FontSize',15)
ylabel('expression sensitivity [a.u.]','FontSize',15)
box off
axis([0.5 2.5 -0.004 inf])


s=subplot(1,2,2);
h=plotSpread(matrix_nt);
set(h{1},'MarkerSize',15,'Color','k')
hold on
boxplot(matrix_nt)
h=findobj(gca,'tag','Outliers'); % Get handles for outlier lines.
set(h,'Marker','none')
h=refline(0,0);
set(h,'LineStyle','--','Color','k')

%calculate AUCs and pvals
[~,~,~,AUC] = perfcurve(ess_oe,NI,true);
p = ranksum(matrix_nt(:,1),matrix_nt(:,2),'tail','left');
text(.8,0.028,['AUC = ' num2str(AUC,2) ', p = ' num2str(p,2)],'FontSize',14)

set(s,'XTick',[1 2],'XTickLabel',{['insensitive, n=' num2str(sum(~ess_oe))],['essential | OE sens., n=' num2str(sum(ess_oe))]},'XTickLabelRotation',45,...
    'YTick',0:.025:0.07,'FontSize',15)
ylabel('noise intolerance [a.u.]','FontSize',15)
box off
axis([0.5 2.5 -0.0075 inf])
print(f,'-depsc2',['results/ESvsNI_essentialgenes.eps'])


f=figure;
s=subplot(1,2,1);
h=plotSpread(matrix_pt1);
set(h{1},'MarkerSize',15,'Color','k')
hold on
boxplot(matrix_pt1)
h=findobj(gca,'tag','Outliers'); % Get handles for outlier lines.
set(h,'Marker','none')

%calculate AUCs and pvals
[~,~,~,AUC] = perfcurve(ess,principal_landscapes_genescores(:,1),true);
p = ranksum(matrix_pt1(:,1),matrix_pt1(:,2),'tail','left');
text(1,max(max(matrix_pt1)),['AUC = ' num2str(AUC,2) ', p = ' num2str(p,2)],'FontSize',14)

set(s,'XTick',[1 2],'XTickLabel',{['non-essential, n=' num2str(sum(~ess))],['essential, n=' num2str(sum(ess))]},'XTickLabelRotation',45,...
    'YTick',0:1:3,'FontSize',15)
ylabel('principal topology 1 loading','FontSize',15)
box off

s=subplot(1,2,2);
h=plotSpread(matrix_pt2);
set(h{1},'MarkerSize',15,'Color','k')
hold on
boxplot(matrix_pt2)
h=findobj(gca,'tag','Outliers'); % Get handles for outlier lines.
set(h,'Marker','none')

%calculate AUCs and pvals
[~,~,~,AUC] = perfcurve(oe,principal_landscapes_genescores(:,2),true);
p = ranksum(matrix_pt2(:,1),matrix_pt2(:,2),'tail','left');
text(1,max(max(matrix_pt2)),['AUC = ' num2str(AUC,2) ', p = ' num2str(p,2)],'FontSize',14)

set(s,'XTick',[1 2],'XTickLabel',{['OE insens., n=' num2str(sum(~oe))],['OE sens., n=' num2str(sum(oe))]},'XTickLabelRotation',45,...
    'YTick',-2:1:1,'FontSize',15)
ylabel('principal topology 2 loading','FontSize',15)
box off

print(f,'-depsc2',['results/PT1_PT2vsDosagesensitivegenes.eps'])

%% compare expression-sensitivity/noise-tolerance against endogenous noise levels

idx_i = ismember(endogenous_validation.Properties.RowNames,genes.ID(genes.GLU_wtExpr > 3 & genes.GLU_wtExpr < 5));

endogenous_validation_limited = endogenous_validation(idx_i,:);


%%% this has to be adapted
noise_idx = 5:7
[r1,p1] = corr(NI',endogenous_validation_limited{:,noise_idx},'rows','pairwise','type','Spearman','tail','left');
[r2,p2] = corr(FNM_pcorr(:,2),endogenous_validation_limited{:,noise_idx},'rows','pairwise','type','Spearman','tail','left');
[r3,p3] = corr(ES',endogenous_validation_limited{:,noise_idx},'rows','pairwise','type','Spearman','tail','left');

Keren2016_dosageSensitivity(Keren2016_dosageSensitivity==Inf) = 10;
[r4,p4] = corr(-Keren2016_dosageSensitivity,endogenous_validation_limited{:,noise_idx},'rows','pairwise','type','Spearman','tail','left');
[r5,p5] = corr(sum(endogenous_validation_limited{:,[3 4]},2)>0,endogenous_validation_limited{:,noise_idx},'rows','pairwise','type','Spearman','tail','left');

r = [r1 ;r2; r3; r4; r5];
p = [p1; p2; p3; p4; p5];

% aggregate pvals using Fisher's method
chi_vals = -2.*log(p);
group_pval = 1 - chi2cdf(sum(chi_vals,2),2*size(p,2))

f=figure;
s=subplot(3,2,[1 3]);
bar(r)
hold on
for i=size(p,2)
    for j=size(p,1)
        if p(j,i) < 0.05
            text(i+(j-2.5)*.17,r(j,i)-0.05,'*','HorizontalAlignment','center')
        end
    end
end
N = sum(~isnan(endogenous_validation_limited{:,[5 6 7]}));
LN = {'NI landscapes','NI raw data','ES','curve Keren','genetic screens'};
set(s,'XTick',1:5,'XTickLabel',LN,'XTickLabelRotation',45,'YTick',-1:.25:1)
ylabel('Spearman rank correlation')
box off


% s=subplot(2,2,2);
% Y = endogenous_validation_limited1{:,7};
% Y(Y==max(Y)) = .4;
% scatter(NI,Y,50,'k','filled')
% hold on
% lsline
% h=refline(0,0);
% set(h,'LineStyle','--','Color','k')
% xlabel('noise intolerance')
% ylabel('total noise diploid')
% set(s,'XTick',-.000:.01:0.045,...
%     'YTick',[-.2:.2:.2 .4],'YTickLabel',[-.2:.2:.2 1.54])
% text(.01,.35,['\rho = ' num2str(r1(3),2) ', p = ' num2str(p1(3),2)])

print(f,'-depsc2',['results/ESvsNI_endogenousnoise.eps'])