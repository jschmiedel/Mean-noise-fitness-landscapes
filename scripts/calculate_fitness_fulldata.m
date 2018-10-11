function calculate_fitness_fulldata()
%calculate fitness from full dataset, including wild-type data

% The doubling time of the wt was 120 minutes in glucose. 
% The read counts for the wt were:
% Glu: time0 - 1421  , time1 - 2253 , time2 - 2679

wt_reads = [1421 2253 2679];
wt_dt = 120/60;

%load genes
load data/genes.mat
%load promoters
load data/promoters.mat
%load reads
GEO_reads0 = readtable('data/Keren2016_GSM2235555_CompetitionDataGlucose.txt','Delimiter','\t');

%get rid of first two cols
GEO_reads1 = GEO_reads0(:,3:end);

%reorder promoters (rows)
[a,b] = ismember(GEO_reads0.PromoterID_SharonEtAl__2012_,promoters.ID);
GEO_reads1(b(a),:) = GEO_reads1(a,:);

varlabels = GEO_reads1.Properties.VariableNames(1:end);
timepoints = cell2mat(cellfun(@(x) str2double(x(end-2:end-1)),varlabels,'UniformOutput',false)');
genes0 = cellfun(@(x) x(1:end-4),varlabels,'UniformOutput',false);

timepoints_unique = unique(timepoints);

%gene lists are preserved in order
for i=1:3
    helper_r = GEO_reads1(:,timepoints==timepoints_unique(i));
    helper_g = genes0(timepoints==timepoints_unique(i));
    [a,b] = ismember(helper_g,genes.ID);
    GEO_reads(:,b(a),i) = helper_r{:,a};
end

%calculate fitness between timepoints
pc=0.1;
tx = 2;
f12 = log2((repmat(wt_reads(1),size(GEO_reads,1),size(GEO_reads,2))./(GEO_reads(:,:,1)+pc)).*...
    ((GEO_reads(:,:,tx)+pc)./repmat(wt_reads(tx),size(GEO_reads,1),size(GEO_reads,2)))) ...
    /((timepoints_unique(tx)-timepoints_unique(1))/wt_dt) + 1;
tx=3;
f13 = log2((repmat(wt_reads(1),size(GEO_reads,1),size(GEO_reads,2))./(GEO_reads(:,:,1)+pc)).*...
    ((GEO_reads(:,:,tx)+pc)./repmat(wt_reads(tx),size(GEO_reads,1),size(GEO_reads,2)))) ...
    /((timepoints_unique(tx)-timepoints_unique(1))/wt_dt) + 1;

%fit linear model
Nisnan = ~isnan(f12) & ~isnan(f13);
P=fit(f12(Nisnan),f13(Nisnan),'poly1');
%and correct one timepoint with it
f12_corr = f12*P.p1+P.p2;


%calculate error per timepoint
pc=0.1; %pseudocount
e12 = sqrt(1./(GEO_reads(:,:,1)+pc) + 1./(GEO_reads(:,:,2)+pc) + 1/wt_reads(1) + 1/wt_reads(2) + 1/400); %add 1/400 term to simulate replicate error of 5%
e13 = sqrt(1./(GEO_reads(:,:,1)+pc) + 1./(GEO_reads(:,:,3)+pc) + 1/wt_reads(1) + 1/wt_reads(3) + 1/400);

%calculate averaged fitness given error of each fitness estimate
f_avg = (f12_corr./(e12.^2)+f13./(e13.^2))./(e12.^-2+e13.^-2);
e123 = sqrt(1./(1./(e12.^2)+1./(e13.^2)));

%plot comparisons between fitness and error from different selection rounds
f=figure;
suptitle('GLU')
subplot(2,2,1)
plot(f12(:),f13(:),'.')
hold on
plot([-1 2],[-1 2]*P.p1 + P.p2,'r:')
xlabel('f12')
ylabel('f13')
axis tight
box off

subplot(2,2,2)
plot(f12_corr(:),f13(:),'.')
hold on
plot([0 1.5],[0 1.5],'r:')
xlabel('f12 corrected')
ylabel('f13')
axis tight
box off

subplot(2,2,3)
loglog(e12(:),e13(:),'.')
xlabel('error f12')
ylabel('error f13')
axis tight
box off

subplot(2,2,4)
semilogy(f_avg(:),e123(:),'.')
xlabel('fitness averaged')
ylabel('error 123')
axis tight
box off

print(f,'-depsc2','results/fitenss_error_calculation.eps')

%calculate fitness, add pseudoread
fitness = array2table(f_avg);
fitness.Properties.VariableNames = genes.ID;
fitness.Properties.RowNames = cellfun(@num2str,num2cell(promoters.ID),'Uni',0);
fitness.Properties.Description = 'GLU';

%calculate error given the read counts
error_fitness = array2table(e123);
error_fitness.Properties.VariableNames = genes.ID;
error_fitness.Properties.RowNames = cellfun(@num2str,num2cell(promoters.ID),'Uni',0);
error_fitness.Properties.Description = 'GLU';

save('data/fitness_t123_glucose.mat','fitness','error_fitness','GEO_reads')