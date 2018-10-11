% calculate surface fits of fitness in mean:noise space
function calculate_surface_fits()

load data/genes.mat
load data/promoters.mat
load data/exclusion_gene_promoters.mat
load data/promoter_expression_GLU_Sharon2014.mat
load data/fitness_t123_glucose.mat

%fitness
Fitness = fitness{:,:};
%only use genes between 3 and 5 mean expression
show_genes = genes.ID(genes.GLU_wtExpr > 3 & genes.GLU_wtExpr < 5);

%set parameters
MODE.Xtra_range = 1;
MODE.weight_xy = 1;
%local grid to evaluate fitness landscape on around each evaluation point
MODE.local_grid_n = 21;

%% estimate dist_weight_xy; the optimal gaussian kernel geometry to minimize RMSE

gridsize = 11;
% pre-define width values of MVN to be assayed
dw_x = logspace(-2,2,gridsize);
dw_y = logspace(-2,2,gridsize);
[hx,hy] = meshgrid(dw_x,dw_y);
dw_xy = [hx(:) hy(:)];

clear RMSE
for J=1:6
    for g=1:numel(show_genes)
        gene = find(strcmp(genes.ID,show_genes{g}));
        Nisnan = promoter_expression.mean>2 & promoter_expression.mean<6 & exclusion_gene_promoters(:,gene)==0;
        
        x=promoter_expression.mean(Nisnan);
        y=promoter_expression.noise(Nisnan);
        z=Fitness(Nisnan,g);
        err_x = promoter_expression.error_mean_modeled(Nisnan);
        err_y = promoter_expression.error_noise_modeled(Nisnan);
        err_z = error_fitness{Nisnan,g};
        
        idx = 1:sum(Nisnan);
        
        %evaluate different distance weights
        for dw = 1:size(dw_xy,1)
            sprintf('%d: %0.1f',J,100*(((g-1)*size(dw_xy,1) + dw)/(numel(show_genes)*size(dw_xy,1))))
            %with dist_weights
            MODE.dist_weight_xy = dw_xy(dw,:);
            
            %10Xval
            RMSE_helper = [];
            for r=0:9
                %define test and training sets by modulus
                test_idx = mod(idx,10) == r;
                training_idx = mod(idx,10) ~= r;
                %evaluate smooth fitness landscape at test promoter locations
                MODE.gridXY = [x(test_idx), y(test_idx)];
                %using as input all training promoters
                F_pred0 = fitness_noise_mean_smoothing_v5(x(training_idx),y(training_idx),z(training_idx),...
                    err_x(training_idx),err_y(training_idx),err_z(training_idx),MODE);
                
                RMSE_helper = [RMSE_helper; ((F_pred0-z(test_idx))./err_z(test_idx)).^2];
            end
            
            RMSE(g,dw,J) = sqrt(mean(RMSE_helper));
            
        end
    end
    
    %update dw_xy grid
    %calculate relative RMSE of each gene in each dw_xy dublet
    %average per dw_xy doublet over all genes
    %pick dw_xy dublet with lowest RMSE
    [~,idx_min_dw_xy] = min(mean(RMSE(:,:,J)./repmat(mean(RMSE(:,:,J),2),1,size(RMSE,2)),1));
    
    min_dw_xy(J+1,:) = log10(dw_xy(idx_min_dw_xy,:));
    
    %update dw_xy
    dw_x(J+1,:) = logspace(min_dw_xy(J+1,1) - (log10(dw_x(J,end))-log10(dw_x(J,1)))/((gridsize-1)*2/3),...
        min_dw_xy(J+1,1) + (log10(dw_x(J,end))-log10(dw_x(J,1)))/((gridsize-1)*2/3),gridsize);
    dw_y(J+1,:) = logspace(min_dw_xy(J+1,2) - (log10(dw_y(J,end))-log10(dw_y(J,1)))/((gridsize-1)*2/3),...
        min_dw_xy(J+1,2) + (log10(dw_y(J,end))-log10(dw_y(J,1)))/((gridsize-1)*2/3),gridsize);
    [hx,hy] = meshgrid(dw_x(J+1,:),dw_y(J+1,:));
    dw_xy = [hx(:) hy(:)];
    
end

MODE.dist_weight_xy = 10.^min_dw_xy(end,:)

figure
for J=1:6

asdf = reshape(mean(RMSE(:,:,J)./repmat(mean(RMSE(:,:,J),2),1,size(RMSE,2)),1),size(dw_y,2),size(dw_y,2));

subplot(3,2,J)
surf(log10(dw_x(J,:)),log10(dw_y(J,:)),asdf) ,
% colorbar
view([-20 75])
xlabel('mean')
ylabel('noise')
zlabel('fitness')
end

%% smooth landscapes with gaussian blurring and numerically estimate gradients

Nisnan = promoter_expression.mean>2 & promoter_expression.mean<6;
MODE.xv = round((min(promoter_expression.mean(Nisnan))-MODE.Xtra_range)*10)/10 : .05 : round((max(promoter_expression.mean(Nisnan))+MODE.Xtra_range)*10)/10;
MODE.n = numel(MODE.xv);
MODE.yv = round((min(promoter_expression.noise(Nisnan))-MODE.Xtra_range)*10)/10 : .025 : round((max(promoter_expression.noise(Nisnan))+MODE.Xtra_range)*10)/10;

MODE = rmfield(MODE,'gridXY');

for gene_idx=1:numel(genes.ID)
    gene_idx
    
    %%%%calculate smoothed landscapes
    
    %exclude data points
    Nisnan = exclusion_gene_promoters(:,gene_idx)==0 & ~isnan(promoter_expression.mean) & promoter_expression.mean>2 & promoter_expression.mean<6; 

    %define data points
    % as averages
    x=promoter_expression.mean(Nisnan);
    y=promoter_expression.noise(Nisnan);
    z=Fitness(Nisnan,gene_idx);
    err_x = promoter_expression.error_mean_modeled(Nisnan);
    err_y = promoter_expression.error_noise_modeled(Nisnan);
    err_z = error_fitness{Nisnan,gene_idx};
    
    smooth_fun.(genes.ID{gene_idx}).F=fitness_noise_mean_smoothing_v5(x,y,z,err_x,err_y,err_z,MODE);
    
    smooth_fun.(genes.ID{gene_idx}).xv = MODE.xv;
    smooth_fun.(genes.ID{gene_idx}).yv = MODE.yv;

    %calculate derivatives
    [smooth_fun.(genes.ID{gene_idx}).Fx,smooth_fun.(genes.ID{gene_idx}).Fy] = ...
        gradient(smooth_fun.(genes.ID{gene_idx}).F,...
        diff(smooth_fun.(genes.ID{gene_idx}).xv([1 2])),diff(smooth_fun.(genes.ID{gene_idx}).yv([1 2])));
    
    smooth_fun.(genes.ID{gene_idx}).Fxx = ...
        gradient(smooth_fun.(genes.ID{gene_idx}).Fx,...
        diff(smooth_fun.(genes.ID{gene_idx}).xv([1 2])),diff(smooth_fun.(genes.ID{gene_idx}).yv([1 2])));
    
end

% correct TUB2
[~,wt_idx] = min(abs(genes.GLU_wtExpr(strcmp(genes.ID,'TUB2'))-smooth_fun.TUB2.xv));
smooth_fun.TUB2.F = smooth_fun.TUB2.F - smooth_fun.TUB2.F(smooth_fun.TUB2.yv==-3,wt_idx) + 1;


%%%%% extract specfic window around each gene's wild-type expression and
%%%%% arrange in one array

gi = fieldnames(smooth_fun);
show_genes = intersect(gi,genes.ID(genes.GLU_wtExpr > 3 & genes.GLU_wtExpr < 5));

%noise levels to evaluate
Neval0 = [-3 -2 -1];
Neval = zeros(numel(Neval0),1);
for i=1:numel(Neval)
    [~,Neval(i)] = min(abs(MODE.yv - Neval0(i)));
end

x_stepsize = smooth_fun.(genes.ID{1}).xv(2) - smooth_fun.(genes.ID{1}).xv(1);
X = repmat(linspace(-1.5,1.5,round(3/x_stepsize)+1)',1,numel(show_genes));
landscape_windowed.X = X;
landscape_windowed.F = NaN([size(X) numel(Neval)]);
landscape_windowed.Fx = NaN([size(X) numel(Neval)]);
landscape_windowed.Fxx = NaN([size(X) numel(Neval)]);
landscape_windowed.Fy = NaN([size(X) numel(Neval)]);
landscape_windowed.F2d = NaN([numel(MODE.yv) size(X,1) size(X,2)]);
landscape_windowed.Fx2d = NaN([numel(MODE.yv) size(X,1) size(X,2)]);
landscape_windowed.Fxx2d = NaN([numel(MODE.yv) size(X,1) size(X,2)]);
landscape_windowed.Fy2d = NaN([numel(MODE.yv) size(X,1) size(X,2)]);

%extract fitness data and derivates at specific noise levels for all genes
%and paste them into array
for   i=1:numel(show_genes)
    gene_idx=find(ismember(genes.ID,show_genes{i}));
    
    [~,wt_idx] = min(abs(smooth_fun.(genes.ID{gene_idx}).xv-genes.GLU_wtExpr(gene_idx)));
    match_gene = max([wt_idx - (size(X,1)-1)/2 1]) : min([wt_idx+(size(X,1)-1)/2 numel(smooth_fun.(genes.ID{gene_idx}).xv)]);
    if numel(match_gene) < size(X,1)
        if wt_idx - (size(X,1)-1)/2 < 1
            match_range = 1:numel(match_gene);
        else
            match_range = (size(X,1)-numel(match_gene)):size(X,1);
        end
    else
        match_range = 1:size(X,1);
    end
    
    for j=1:numel(Neval)
        landscape_windowed.F2d(:,match_range,i) = smooth_fun.(genes.ID{gene_idx}).F(:,match_gene);
        landscape_windowed.Fx2d(:,match_range,i) = smooth_fun.(genes.ID{gene_idx}).Fx(:,match_gene);
        landscape_windowed.Fxx2d(:,match_range,i) = smooth_fun.(genes.ID{gene_idx}).Fxx(:,match_gene);
        landscape_windowed.Fy2d(:,match_range,i) = smooth_fun.(genes.ID{gene_idx}).Fy(:,match_gene);
        
        landscape_windowed.F(match_range,i,j) = smooth_fun.(genes.ID{gene_idx}).F(Neval(j),match_gene);
        landscape_windowed.Fx(match_range,i,j) = smooth_fun.(genes.ID{gene_idx}).Fx(Neval(j),match_gene);
        landscape_windowed.Fxx(match_range,i,j) = smooth_fun.(genes.ID{gene_idx}).Fxx(Neval(j),match_gene);
        landscape_windowed.Fy(match_range,i,j) = smooth_fun.(genes.ID{gene_idx}).Fy(Neval(j),match_gene);
    end
    
end

save('data/FNM_smooth_fun_grid05_025.mat','smooth_fun','landscape_windowed','Nisnan','MODE','Neval')



