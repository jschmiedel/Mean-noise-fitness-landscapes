function fitness_landscape_PCA()

load data/genes.mat
load data/promoters.mat
load data/exclusion_gene_promoters.mat

load data/FNM_smooth_fun_grid05_025.mat
show_genes = genes.ID(genes.GLU_wtExpr > 3 & genes.GLU_wtExpr < 5);

%normalize fitness by fitness at wild-type mean and low noise
F2d_norm = landscape_windowed.F2d(:,:,:)-...
    repmat(landscape_windowed.F2d(Neval(1),(size(landscape_windowed.F2d,2)-1)/2,:),size(landscape_windowed.F2d,1),size(landscape_windowed.F2d,2));

%take only data within noise=[-3,-1]
rel_fit = (squeeze(reshape(F2d_norm(MODE.yv>=-3 & MODE.yv<=-1,:,:),sum(MODE.yv>=-3 & MODE.yv<=-1)*size(landscape_windowed.F2d,2),1,size(landscape_windowed.F2d,3))))';


%% PCA
[COEFF0, SCORE0, LATENT] = pca(rel_fit,'Centered',true);
COEFF = -COEFF0;
SCORE = -SCORE0;

%estimate PC coefficients for mean landscape
fun2 = @(y) sqrt(sum((COEFF(:,[1 2])*y - mean(rel_fit)').^2));
y0 = zeros(2,1);
PCzero = fminsearch(fun2,y0);

principal_landscapes = reshape(COEFF,81,61,32);
principal_landscapes_genescores = SCORE(:,1:2) + repmat(PCzero',size(SCORE,1),1); %correct loadings for mean landscape loading
principal_landscapes_varexp = LATENT;

save(['data/FNM_smooth_fun_grid05_025_landscape_PCA.mat'],'principal_landscapes','principal_landscapes_genescores','principal_landscapes_varexp','PCzero')

%%
f=figure;
subplot(3,2,1)
plot(principal_landscapes_genescores(:,1),principal_landscapes_genescores(:,2),'.')
for i=1:numel(show_genes)
    text(principal_landscapes_genescores(i,1),principal_landscapes_genescores(i,2),show_genes(i))
end
hline(0)
refline(1,0)
refline(-1,0)
box off
axis([0 ...
    max(principal_landscapes_genescores(:,1))+0.1*range(principal_landscapes_genescores(:,1)) ...
    min(principal_landscapes_genescores(:,2))-0.05*range(principal_landscapes_genescores(:,2)) ...
    max(principal_landscapes_genescores(:,2))+0.1*range(principal_landscapes_genescores(:,2))])
xlabel(['PC 1 loading (' num2str(round(LATENT(1)/sum(LATENT)*100,1)) '% variance ex.)'])
ylabel(['PC 2 loading (' num2str(round(LATENT(2)/sum(LATENT)*100,1)) '% variance ex.)'])


s=subplot(3,3,3);
bar(LATENT/sum(LATENT))
set(s,'XTick',[1 5:5:30],'YTick',[0 0.5 1])
axis([0 33 0 1])
xlabel('principal component')
ylabel('% variance explained')


s=subplot(3,3,4);
p=pcolor(-1.5:0.05:1.5,MODE.yv(MODE.yv>=-3 & MODE.yv<=-1),principal_landscapes(:,:,1));
set(p, 'EdgeColor', 'none');
hold on
V = -.055:.005:.025; 
[~,c]=contour(-1.5:0.05:1.5,MODE.yv(MODE.yv>=-3 & MODE.yv<=-1),principal_landscapes(:,:,1),V);
set(c,'LineColor','black');

% colorbar
caxis(V([1 end]))
set(s,'XTick',[-1 0 1],'YTick',[-3 -2 -1])
colormap('jet')
title('PC 1')
xlabel('mean expression')
ylabel('noise')

%%%%plot to create colorbar
s=subplot(3,3,6);
contourf(-1.5:0.05:1.5,MODE.yv(MODE.yv>=-3 & MODE.yv<=-1),principal_landscapes(:,:,1),V);
caxis(V([1 end]))
colorbar
title('ignore this; just for colbar')
%%%%


s=subplot(3,3,5);
p=pcolor(-1.5:0.05:1.5,MODE.yv(MODE.yv>=-3 & MODE.yv<=-1),principal_landscapes(:,:,2));
set(p, 'EdgeColor', 'none');
hold on
[~,c]=contour(-1.5:0.05:1.5,MODE.yv(MODE.yv>=-3 & MODE.yv<=-1),principal_landscapes(:,:,2),V);
set(c,'LineColor','black');
caxis(V([1 end]))
set(s,'XTick',[-1 0 1],'YTick',[-3 -2 -1])
title('PC 2')
xlabel('mean expression')
ylabel('noise')


s=subplot(3,3,7);
p=pcolor(-1.5:0.05:1.5,MODE.yv(MODE.yv>=-3 & MODE.yv<=-1),principal_landscapes(:,:,1) + principal_landscapes(:,:,2));
set(p, 'EdgeColor', 'none');
hold on
[~,c]=contour(-1.5:0.05:1.5,MODE.yv(MODE.yv>=-3 & MODE.yv<=-1),principal_landscapes(:,:,1) + principal_landscapes(:,:,2),V);
set(c,'LineColor','black');
caxis(V([1 end]))
set(s,'XTick',[-1 0 1],'YTick',[-3 -2 -1])
title('PC 1 + PC 2')
xlabel('mean expression')
ylabel('noise')

s=subplot(3,3,8);
p=pcolor(-1.5:0.05:1.5,MODE.yv(MODE.yv>=-3 & MODE.yv<=-1),principal_landscapes(:,:,1) - principal_landscapes(:,:,2));
set(p, 'EdgeColor', 'none');
hold on
[~,c]=contour(-1.5:0.05:1.5,MODE.yv(MODE.yv>=-3 & MODE.yv<=-1),principal_landscapes(:,:,1) - principal_landscapes(:,:,2),V);
set(c,'LineColor','black');
caxis(V([1 end]))
set(s,'XTick',[-1 0 1],'YTick',[-3 -2 -1])
title('PC 1 - PC 2')
xlabel('mean expression')
ylabel('noise')

s=subplot(3,3,9);
F_mean = reshape(mean(rel_fit),81,61);
F_mean_reconstruct = PCzero(1)*principal_landscapes(:,:,1) + PCzero(2)*principal_landscapes(:,:,2);

p=pcolor(-1.5:0.05:1.5,MODE.yv(MODE.yv>=-3 & MODE.yv<=-1),F_mean-F_mean_reconstruct);
set(p, 'EdgeColor', 'none');
hold on
[~,c]=contour(-1.5:0.05:1.5,MODE.yv(MODE.yv>=-3 & MODE.yv<=-1),F_mean-F_mean_reconstruct,V);
set(c,'LineColor','black');
caxis(V([1 end]))
set(s,'XTick',[-1 0 1],'YTick',[-3 -2 -1])
title(['diff mean vs. mean-reconst, R^2 = ' num2str(corr(F_mean_reconstruct(:),F_mean(:))^2,2)])
xlabel('mean expression')
ylabel('noise')

print(f,'-depsc2','-painters','-loose','results/landscape_PCA.eps')