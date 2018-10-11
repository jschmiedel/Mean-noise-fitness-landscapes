function plot_FNM_scatter_rawdata_v2()
%% plot fitness:noise:mean 'raw' data for (example) genes

load data/genes.mat
load data/promoters.mat
load data/exclusion_gene_promoters.mat

load data/promoter_expression_GLU_Sharon2014.mat
load data/fitness_t123_glucose.mat

load data/FNM_smooth_fun_grid05_025.mat
load data/FNM_smooth_fun_grid05_025_landscape_PCA.mat

Fitness = fitness{:,:};

show_genes = genes.ID(genes.GLU_wtExpr > 3 & genes.GLU_wtExpr < 5);
%% first ALL genes

%colorscale
colorscale = 0.8:.01:1.1;

f=figure;
pidx = 1;

%plot colorbar
subplot(3,3,1);
gene = find(strcmp(genes.ID,show_genes{1}));
Nisnan = exclusion_gene_promoters(:,gene)==0;
scatter(promoter_expression.mean(Nisnan),promoter_expression.noise(Nisnan),sqrt(1./error_fitness{Nisnan,gene})*10,fitness{Nisnan,gene},'filled')
colormap('jet')
colorbar
caxis([colorscale(1) colorscale(end)])

%plot all genes
for i=1:numel(show_genes)
    
    if i+1 > pidx*9
        if pidx == 1
            print(f,'-dpsc2','results/raw_data_FNM_allgenes.ps')
        else
            print(f,'-dpsc2','results/raw_data_FNM_allgenes.ps','-append')
        end
        %         close(f)
        pidx = pidx+1;
        f=figure;
    end
    
    
    
    s=subplot(3,3,i+1-(pidx-1)*9);
    
    gene = find(strcmp(genes.ID,show_genes{i}));
    if strcmp(show_genes{i},'TUB2') %correct TUB2
        Fitness(Nisnan,gene) = fitness{Nisnan,gene} + 0.1019;
    end
    
    Nisnan = exclusion_gene_promoters(:,gene)==0;
    errbar(promoter_expression.mean(Nisnan),promoter_expression.noise(Nisnan),promoter_expression.error_noise(Nisnan),'Color',ones(3,1)*.5)
    hold on
    errbar(promoter_expression.mean(Nisnan),promoter_expression.noise(Nisnan),promoter_expression.error_mean(Nisnan),'horiz','Color',ones(3,1)*.5)
    scatter(promoter_expression.mean(Nisnan),promoter_expression.noise(Nisnan),sqrt(1./error_fitness{Nisnan,gene})*10,Fitness(Nisnan,gene),'filled')
    plot(ones(2,1)*genes.GLU_wtExpr(gene),[-3.5 -0.5],'k:','LineWidth',0.5)
    plot(ones(2,1)*(genes.GLU_wtExpr(gene)-1.5),[-3.5 -0.5],'r:','LineWidth',0.5)
    plot(ones(2,1)*(genes.GLU_wtExpr(gene)+1.5),[-3.5 -0.5],'r:','LineWidth',0.5)
    colormap('jet')
    caxis([colorscale(1) colorscale(end)])
    title(show_genes{i},'FontSize',10)
    set(s,'XTick',2:6,'YTick',-3:-1,'FontSize',9)
    xlabel('expression (log2) [a.u.]','FontSize',8)
    ylabel('noise [log2(CV)]','FontSize',8)
    axis([1.5 6.5 -3.5 -0.5])
    box off
end

print(f,'-dpsc2','results/raw_data_FNM_allgenes.ps','-append')

%% plot 3 example genes

example_genes = {'TOP1','RPN8','TUB2'};

%colorscale
colorscale = 0.8:.01:1.1;

f=figure;
pidx = 1;

%plot colorbar
subplot(3,3,4);
gene = find(strcmp(genes.ID,show_genes{1}));
Nisnan = exclusion_gene_promoters(:,gene)==0;
scatter(promoter_expression.mean(Nisnan),promoter_expression.noise(Nisnan),sqrt(1./error_fitness{Nisnan,gene})*10,fitness{Nisnan,gene},'filled')
colormap('jet')
colorbar
caxis([colorscale(1) colorscale(end)])

%plot all genes
for i=1:numel(example_genes)
    
    s=subplot(3,3,i);
    
    gene = find(strcmp(genes.ID,example_genes{i}));
    if strcmp(example_genes{i},'TUB2') %correct TUB2
        Fitness(Nisnan,gene) = fitness{Nisnan,gene} + 0.1019;
    end
    Nisnan = exclusion_gene_promoters(:,gene)==0;
    errbar(promoter_expression.mean(Nisnan),promoter_expression.noise(Nisnan),promoter_expression.error_noise(Nisnan),'Color',ones(3,1)*.5)
    hold on
    errbar(promoter_expression.mean(Nisnan),promoter_expression.noise(Nisnan),promoter_expression.error_mean(Nisnan),'horiz','Color',ones(3,1)*.5)
    scatter(promoter_expression.mean(Nisnan),promoter_expression.noise(Nisnan),sqrt(1./error_fitness{Nisnan,gene})*10,Fitness(Nisnan,gene),'filled')
    plot(ones(2,1)*genes.GLU_wtExpr(gene),[-3.5 -0.5],'k:','LineWidth',0.5)
    plot(ones(2,1)*(genes.GLU_wtExpr(gene)-1.5),[-3.5 -0.5],'r:','LineWidth',0.5)
    plot(ones(2,1)*(genes.GLU_wtExpr(gene)+1.5),[-3.5 -0.5],'r:','LineWidth',0.5)
    colormap('jet')
    caxis([colorscale(1) colorscale(end)])
    title(example_genes{i},'FontSize',8)
    set(s,'XTick',2:6,'YTick',-3:-1,'FontSize',8)
    xlabel('mean expression (log2) [a.u.]','FontSize',8)
    ylabel('expression noise [log2(CV)]','FontSize',8)
    axis([1.5 6.5 -3.5 -0.5])
    box off
    
    
    s=subplot(3,3,i+6);
    p=pcolor(smooth_fun.(example_genes{i}).xv,smooth_fun.(example_genes{i}).yv,smooth_fun.(example_genes{i}).F);
    set(p, 'EdgeColor', 'none');
    hold on
    [~,c]=contour(smooth_fun.(example_genes{i}).xv,smooth_fun.(example_genes{i}).yv,smooth_fun.(example_genes{i}).F,colorscale);
    set(c,'LineColor','black');
    hold on
    scatter(promoter_expression.mean(Nisnan),promoter_expression.noise(Nisnan),20,.75*ones(sum(Nisnan),3),'filled')
    plot(ones(2,1)*genes.GLU_wtExpr(gene),[-3.5 -0.5],'k:','LineWidth',0.5)
    plot(ones(2,1)*(genes.GLU_wtExpr(gene)+1.5),[-3.5 -0.5],'r:','LineWidth',0.5)
    plot(ones(2,1)*(genes.GLU_wtExpr(gene)-1.5),[-3.5 -0.5],'r:','LineWidth',0.5)
    caxis([colorscale(1) colorscale(end)])
    % colorbar
    colormap('jet')
    axis([2 6 -3 -1])
    set(s,'XTick',2:6,'YTick',-3:-1,'FontSize',8)
    xlabel('mean expression (log2) [a.u.]','FontSize',8)
    ylabel('expression noise [log2(CV)]','FontSize',8)
    box off
    
end

print(f,'-depsc2','results/raw_data_FNM_examplegenes.eps')

%% show for one example gene the transformation from raw FNM data to smoothed landscape

GENE = 'TUB2';
gene = find(strcmp(genes.ID,GENE));


f=figure;
s=subplot(2,2,1);
x=linspace(1.5,6.5,80);
y=linspace(-3.5,-0.5,80);
[X1,Y1] = meshgrid(x,y);
D = sqrt((X1 - 4).^2 + (Y1 + 1.5).^2);
% W = exp(-(D.^2));
W = mvnpdf([X1(:) Y1(:)],[4 -2],MODE.dist_weight_xy);
W = reshape(W,numel(y),numel(x));
W = W/max(max(W));

h=pcolor(x,y,W);
set(h,'EdgeColor','none')
% imagesc(x,y,W)
hold on
scatter(promoter_expression.mean(Nisnan),promoter_expression.noise(Nisnan),20,.75*ones(sum(Nisnan),3),'filled')
[~,maxxi] = max(W);
[~,maxyi] = max(W');
plot(4,-2,'k.','MarkerSize',20);
[c,h]=contour(x,y,W,[.25 .5 .75 .9],'k');
clabel(c,h,'FontSize',10)
caxis([0 1])
white = 225;
Cmap = ([linspace(0,white,31) white*ones(1,2) linspace(white,237,31);...
    linspace(114,white,31) white*ones(1,2) linspace(white,177,31);...
    linspace(189,white,31) white*ones(1,2) linspace(white,32,31)]'/255);
Cmap = ([linspace(255,0,301);...
    linspace(255,114,301);...
    linspace(255,189,301)]'/255);
colormap(Cmap)
colorbar
axis([2 6 -3 -1])
title('gaussian smoothing','FontSize',16)
set(s,'XTick',2:6,'YTick',-3:-1,'FontSize',12,'YDir','normal')
xlabel('expression (log2) [a.u.]','FontSize',14)
ylabel('noise [log2(CV)]','FontSize',14)
box off
print(f,'-depsc2','-painters','-loose',['results/setup_smoothfun.eps'])

%% expression fitness function in comparison to mean-noise fitness landscape
GENE = 'TUB2';
gene = find(strcmp(genes.ID,GENE));

if strcmp(GENE,'TUB2') %correct TUB2
    Fitness(Nisnan,gene) = fitness{Nisnan,gene} + 0.1019;
end
Nisnan = exclusion_gene_promoters(:,gene)==0;

f=figure;
s=subplot(3,3,[1 2]);
scatter(promoter_expression.mean(Nisnan),Fitness(Nisnan,gene),35,'k','filled')
hold on
r=refline(0,1);
plot(ones(2,1)*genes.GLU_wtExpr(gene),[0.5 1.125],'r--','LineWidth',0.5)
set(r,'Color',[0.5 0.5 0.5],'LineWidth',0.5)
axis([2 6 0.5 1.125])
set(s,'Ytick',0.4:0.2:1,'FontSize',14)
xlabel('mean expression [log2]','FontSize',18)
ylabel('fitness','FontSize',18)

s=subplot(3,5,[6 14]);
colorscale = 0.6:.01:1.1;
hold on
scatter(promoter_expression.mean(Nisnan),promoter_expression.noise(Nisnan),50,Fitness(Nisnan,gene),'filled')
plot(ones(2,1)*genes.GLU_wtExpr(gene),[-3.5 -0.5],'r--','LineWidth',0.5)
colormap('jet')
colorbar
caxis([colorscale(1) colorscale(end)])
set(s,'XTick',2:6,'YTick',-3:-1,'FontSize',14)
xlabel('mean expression [log2]','FontSize',18)
ylabel('expression noise [log2]','FontSize',18)
axis([2 6 -3.1 -0.9])
box off
print(f,'-depsc2','-painters','-loose',['results/setup_example_2dvs3d_' GENE '.eps'])

%% plot all landscapes

f=figure;
for i=1:numel(show_genes)
    s = subplot(7,5,i+2);
    F2d = landscape_windowed.F2d(smooth_fun.(show_genes{i}).yv >= -3 & smooth_fun.(show_genes{i}).yv <=-1,:,i);
    F2d = F2d-F2d(1,31);
    p=pcolor(-1.5:0.05:1.5,-3:0.025:-1,F2d);
    set(p, 'EdgeColor', 'none');
    hold on
    V = -0.2:0.01:0.1;
    [~,c]=contour(-1.5:0.05:1.5,-3:0.025:-1,F2d,V);
    set(c,'LineColor','black');
    
    % axis([2 6 -3 -1])
    axis([-1.5 1.5 -3 -1])
    caxis(V([1 end]))
    set(s,'XTick',[],'YTick',[])
    colormap('jet')
    title([show_genes{i},' (',num2str(round(principal_landscapes_genescores(i,1),1)),';',num2str(round(principal_landscapes_genescores(i,2),1)),')'],'FontSize',10)
    xlabel('')
    ylabel('')
end

s=subplot(7,5,1);
axis([-1.5 1.5 -3 -1])
set(s,'Xtick',-1:1:1,'XTickLabel',{'-1','wt','1'},'YTick',-3:1:-1,'FontSize',12)
xlabel('mean','FontSize',12)
ylabel('noise','FontSize',12)

% this was manually printed to PDF after resizing

%% plot ENO2 and RPL3 landscapes (transition in noise-fitness effects)

X{1} = 'ENO2';
X{2} = 'RPL3';

f=figure;
for i=1:numel(X)
    s = subplot(2,2,i);
    F2d = smooth_fun.(X{i}).F(smooth_fun.(X{i}).yv >= -3 & smooth_fun.(X{i}).yv <=-1,...
        smooth_fun.(X{i}).xv >= 2 & smooth_fun.(X{i}).xv <=6);
    
    p=pcolor(2:0.05:6,-3:0.025:-1,F2d);
    set(p, 'EdgeColor', 'none');
    hold on
    V = 0.8:0.01:1.1;
    [~,c]=contour(2:0.05:6,-3:0.025:-1,F2d,V);
    set(c,'LineColor','black');
    
    axis([2 6 -3 -1])
    caxis(V([1 end]))
    set(s,'XTick',2:6,'YTick',-3:-1)
    colormap('jet')
    title([X{i},', wildtype expr. = ',num2str(round(genes.GLU_wtExpr(strcmp(genes.ID,X{i})),1))],'FontSize',12)
    xlabel('mean expression (log2)')
    ylabel('noise (log2)')
end

print(f,'-depsc2','-painters','-loose',['results/examples_noisefitness_transition_ENO2_RPL3.eps'])