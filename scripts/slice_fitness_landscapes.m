function slice_fitness_landscapes()
%% slice cluster-average fitness landscape in mean and noise direction

load data/genes.mat
load data/promoters.mat
load data/exclusion_gene_promoters.mat

load data/promoter_expression_GLU_Sharon2014.mat
load data/fitness_t123_glucose.mat

load data/FNM_smooth_fun_grid05_025.mat
load data/FNM_smooth_fun_grid05_025_landscape_PCA.mat

pc_combos = [1 0;
    0 1;
    1 1;
    1 -1];

pc_combo_labels = {'PC 1','PC 2','PC1 + PC2','PC1 - PC2'};

%% fitness as function of expression at noise=-2 for all clusters

Colu = linspecer(size(pc_combos,1));
f=figure;
for u=1:size(pc_combos,1)
    h = pc_combos(u,1) * principal_landscapes(:,:,1) + pc_combos(u,2) * principal_landscapes(:,:,2);
    x = landscape_windowed.X(:,1);
    l(u) = plot(x(:,1),h((size(principal_landscapes,1)-1)/2*(j-1)+1,:),'-','Color',Colu(u,:),'LineWidth',1);
    hold on
end
r=refline(0,0);
set(r,'Color','k','LineStyle','--','LineWidth',0.5)
box off
xlabel('expression rel. to WT [log2]','FontSize',14)
ylabel('fitness rel. to WT','FontSize',14)
axis([-1.5 1.5 -0.035 0.035])
set(gca,'XTick',[-1 0 1],'YTick',-0.04:.01:0.04,'FontSize',12)
h=legend(l,pc_combo_labels,'Location','Best');
set(h,'FontSize',12)
title(['log2(noise) = ' num2str(MODE.yv(Neval(j)))])

print(f,'-depsc2','-painters','-loose','results/landscape_slice_Fx_allPCs.eps')

%% fitness as function of noise at wild-type expression

f=figure;
Col3 = linspecer(size(pc_combos,1));
x = MODE.yv(MODE.yv >= -3 & MODE.yv <= -1);
clear l
for u=1:size(pc_combos,1)
    hy = pc_combos(u,1) * principal_landscapes(:,:,1) + pc_combos(u,2) * principal_landscapes(:,:,2);
    l(u)=plot(x,hy(:,(size(principal_landscapes,2)-1)/2+1),'-','Color',Col3(u,:),'LineWidth',2);
    hold on
end
r=refline(0,0);
set(r,'Color','k','LineStyle','--')
legend(l,pc_combo_labels,'Location','EastOutside')
box off
set(gca,'XTick',[-3 -2 -1],'YTick',-0.05:.01:0.01,'FontSize',12)
xlabel('noise [log2(CV)]','FontSize',16)
ylabel('rel. fitness','FontSize',16)
title('fitness at WT expr.','Fontsize',14)
print(f,'-depsc2','-painters','-loose','results/landscape_slice_Fvsnoise_wt_allPCs.eps')