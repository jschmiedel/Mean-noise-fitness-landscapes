function evolve_fitness_landscapes()

%% stochastic evolutionary trajectories on principal topologies & peaked
%landscape

load data/FNM_smooth_fun_grid05_025.mat
load data/FNM_smooth_fun_grid05_025_landscape_PCA.mat

%vector of possible mutations; burst size mutations (only changing mean: 1,2) and burst
%frequency mutations (changing both mean and noise: 3,4)
mut_vector = [0 1;
    0 -1; ...
    -1 1; ...
    1 -1];

%mutation probability; all mutations equally likely
mut_p{1} = 1./abs(mut_vector(:,2));

%are burst size mutations less likely than burst frequency mutations?
% size mutations 10 times less likely
mut_p{2} = 1./abs(mut_vector(:,2));
mut_p{2}(1:2) = mut_p{2}(1:2)/10;

% size mutations 10 times more(!) likely
mut_p{3} = 1./abs(mut_vector(:,2));
mut_p{3}(3:4) = mut_p{3}(3:4)/10;

landscape_size = size(principal_landscapes(:,:,1));

%define different population start points in mean-noise space
start_positions = [(2:(landscape_size(1)-1))' repmat(2,landscape_size(1)-2,1); %low mean, diff noise
    (2:(landscape_size(1)-1))' repmat(landscape_size(2)-1,landscape_size(1)-2,1); %high mean, diff noise
    repmat(landscape_size(1)-1,landscape_size(2)-2,1) (2:(landscape_size(2)-1))']; %high noise, diff mean

pc_combos = [1 0;
    0 1;
    1 1.2];

pc_combo_labels = {'PC 1','PC 2','PC1 + PC2'};

for m=1:numel(mut_p)
    for j=1:size(pc_combos,1)
        F = principal_landscapes(:,:,1) * pc_combos(j,1) + principal_landscapes(:,:,2) * pc_combos(j,2);
        for s=1:size(start_positions,1)
            %noise-mean-time vector
            nmt_vector = [start_positions(s,1) start_positions(s,2) 0];
            
            %gillepsie algorithm
            boundary_not_yet_reached = true;
            while boundary_not_yet_reached
                
                helper = repmat(nmt_vector(end,[1 2]),4,1) + mut_vector;
                delta_fitness = F(sub2ind(landscape_size,helper(:,1),helper(:,2))) - F(nmt_vector(end,1),nmt_vector(end,2));
                
                p = mut_p{m} .* arrayfun(@(x) max([0 x]),delta_fitness);
                r1 = rand;
                tau = -log(r1)/sum(p);
                
                nmt_vector(end+1,3) = nmt_vector(end,3)+tau;
                r2 = rand;
                case_r2 = find(cumsum(p)/sum(p)-r2>0,1,'first');
                nmt_vector(end,[1 2]) = nmt_vector(end-1,[1 2]) + mut_vector(case_r2,:);
                
                if nmt_vector(end,1) == 1 || nmt_vector(end,2) == 1 || nmt_vector(end,1) == landscape_size(1) || nmt_vector(end,2) == landscape_size(2)
                    boundary_not_yet_reached = false;
                end
            end
            evolutionary_path{j,s,m} = nmt_vector;
        end
    end
end

save(['results/FNM_smooth_fun_grid05_025_evoLandscape_gillespie_all.mat'],'evolutionary_path','pc_combos','mut_vector','mut_p','start_positions')