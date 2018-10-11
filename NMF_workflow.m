%analyse noise-fitness data; workflow script


%define source folder
FOLDER_prefix = 'where_is_the_github_repository?';

%add script path
path(path,[FOLDER_prefix 'scripts/'])
path(path,[FOLDER_prefix 'scripts/supplemental'])

cd(FOLDER_prefix)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% general data processing %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%% genes
define_genes()

%% %%% promoters
define_promoters()

%% %%% promoter expression: mean & noise
promoter_expression_Sharon2014()

%% %%% calculate fitness
calculate_fitness_fulldata()

%% %%% exclude promoter-gene combinations
exclude_promoter_gene_combos()

%% %%% collect data from large-scale screens on expression-sensitivity and noise of genes

assemble_largescale_data_for_validation()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% analyses %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate surface fits in FNM space
calculate_surface_fits()
%the results from this are given as data/FNM_smooth_fun_grid05_025.mat

%% PCA on surface fits 
fitness_landscape_PCA()

%% slice principal topologies in mean and noise direction (fig2)
slice_fitness_landscapes()

%% calculate expression-sensitivity and noise intolerance on landscapes (fig 3)
calculate_exprsens_noiseintol_landscapes()

%% analyse fitness-sensitivity to mean/noise changes on (cluster-average) fitness landscapes (fig4)

fitness_sensitivity_noise_mean_changes()

%% evolution on fitness landscapes
% 1) simulate evolutionary trajectories
evolve_fitness_landscapes()

%% 2) visualize trajectories on landscapes and their dynamics
evaluate_evolved_fitness_landscapes()

%% 3) epistasis of promoter mutations on fitness landscape
epistasis_on_fitness_landscapes_v2()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% plot %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot raw scatter and transition to smoothed landscapes

%plot raw data, smoothed landscapes and PC reconstruction for all genes
plot_FNM_scatter_rawdata_v2()