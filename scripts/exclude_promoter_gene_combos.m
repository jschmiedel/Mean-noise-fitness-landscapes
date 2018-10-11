function exclude_promoter_gene_combos()

load data/promoters.mat
load data/genes.mat

exclusion_gene_promoters = sparse(numel(promoters.ID),numel(genes.ID));

% check error between mean estimates of Keren2016 and Sharon2014
load('data/promoter_expression_GLU_Sharon2014.mat')
x_mean = abs((promoter_expression.mean-promoters.GLU_expression_Keren2016))>0.5;

% exclude promoter-gene combos with too high error of fitness estimate
load('data/promoter_expression_GLU_Sharon2014.mat')
load('data/fitness_t123_glucose.mat')
x_fitness = nanmedian(error_fitness{:,:},2) > 0.1;

% exclude promoters below 2 or above 6 expression
x_expr = ~(promoter_expression.mean>2 & promoter_expression.mean<6);

exclusion_gene_promoters(x_mean|x_fitness|x_expr,:) = 1;

['promoters with large mean expression differences between the two studies: ' num2str(sum(x_mean))]
['+promoters with large fitness errors: ' num2str(sum(~x_mean & x_fitness))]
['++promoters outside of core expression region: ' num2str(sum(~x_mean & ~x_fitness & x_expr))]

sum(~x_mean & ~x_expr & ~x_fitness)

%% per gene, discard promoters if gene is TF and has binding site in promoter

% how many promoer:gene combos are present?
helper = sum(sum(~exclusion_gene_promoters));


Sharon2012_TFmotives = readtable('data/Sharon2012_TableS1.xlsx'); %% copy-pasted from Sharon2012 supplement
Sharon2012_PromSeqs0 = readtable('data/Sharon2012_TableS3.xlsx'); %% copy-pasted from Sharon2012 supplement

%first clear nonsense lines in Sharon2012_promoters
Sharon2012_TFmotives = Sharon2012_TFmotives(strncmp(Sharon2012_TFmotives.Type,'TF',2),:);
%trim white-spaces and cut gene names
%gene name text bits to remove
remove_bits = {'_BARAK','AllPSSMs_','_Hovring','_HellauerLEU2','_2','Core','_flipped','_v1','_v2','Giniger_','_Site1in2','_Site3in4','_17merin4','_Site1','_Site2','_Site3','_Site4','_17mer'};
for i=1:size(Sharon2012_TFmotives,1)
    Sharon2012_TFmotives.Name{i}(strfind(Sharon2012_TFmotives.Name{i},' ')) = '';
    for j=1:numel(remove_bits)
        match = strfind(Sharon2012_TFmotives.Name{i},remove_bits{j});
        if ~isempty(match)
            Sharon2012_TFmotives.Name{i}(match:match-1+length(remove_bits{j})) = '';
        end
    end
    Sharon2012_TFmotives.Sequence{i}(strfind(Sharon2012_TFmotives.Sequence{i},' ')) = '';
end

%exclude MSN2 and MSN4 from analysis, because their GGGG motif is found in every promoter
Sharon2012_TFmotives = Sharon2012_TFmotives(~strncmp(Sharon2012_TFmotives.Name,'MSN',3),:);

% map TF binding sites to promoter sequences
%get sequences of promoters
[a,~]=ismember(Sharon2012_PromSeqs0.LibraryID,promoters.ID);
Sharon2012_PromSeqs = Sharon2012_PromSeqs0.OligoSequence(a,:);

for i=1:size(Sharon2012_TFmotives,1)
    match_gene = find(strcmp(Sharon2012_TFmotives.Name{i},genes.ID));
    if ~isempty(match_gene)
        for j=1:length(Sharon2012_PromSeqs)
            match_motif = strfind(Sharon2012_PromSeqs{j},Sharon2012_TFmotives.Sequence{i});
            if ~isempty(match_motif)
                exclusion_gene_promoters(j,match_gene) = true;
            end
        end
    end
end
  
['+++gene:promoters combos with potential self-regulation: ' num2str(helper-sum(sum(~exclusion_gene_promoters)))]

%% plot promoter choices
f=figure;

s=subplot(2,2,1);
[h,x]=hist(genes.GLU_wtExpr,0.5:7.5);
bar(x,h)
hold on
plot(ones(2,1)*3,[0 25],'k--')
plot(ones(2,1)*5,[0 25],'k--')
set(s,'XTick',2:1:6,'YTick',0:10:30)
axis([-.5 8 0 25])
ylabel('# genes')
xlabel('wild-type expression (log2) [a.u.]')
box off

s=subplot(2,2,2);
[h,x]=hist(nanmedian(error_fitness{:,:},2),-0.0125:0.025:0.2625);
bar(x,h)
hold on
plot(ones(2,1)*.1,[0 25],'k--')
set(s,'XTick',0:0.05:0.25,'YTick',0:20:60)
axis([0 0.25 0 75])
ylabel('# promoters')
xlabel('median error in fitness measure [a.u.]')
box off


s=subplot(2,2,[3 4]);
G = 4*ones(numel(promoter_expression.mean),1);
G(x_expr) = 3;
G(x_fitness) = 2;
G(x_mean) = 1;

gscatter(promoter_expression.mean,promoter_expression.noise,G,'rbgk','....o',[],0)
hold on
plot(ones(2,1)*2,[-3.3 -0.5],'k--')
plot(ones(2,1)*6,[-3.3 -0.5],'k--')
xlabel('mean expression (log2) [a.u.]')
ylabel('expression noise log2(CV)')
l=legend(['mean diff(Sharon14:Keren16) > 0.5, n=' num2str(sum(G==1))],...
    ['median(fitness error) > .1, n=' num2str(sum(G==2))],...    
    ['expression off, n=' num2str(sum(G==3))],...
    ['used, n=' num2str(sum(G==4))],'Location','EastOutside');
set(s,'XTick',2:6,'YTick',-3:-1)
axis([-.5 8 -3.3 -0.5])

box off

print(f,'-depsc2','results/setup_promoter_use.eps')
%%
save('data/exclusion_gene_promoters.mat','exclusion_gene_promoters')