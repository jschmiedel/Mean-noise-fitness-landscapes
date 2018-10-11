function evaluate_evolutionary_trajectories()

%%
load data/FNM_smooth_fun_grid05_025.mat
load data/FNM_smooth_fun_grid05_025_landscape_PCA.mat
load results/FNM_smooth_fun_grid05_025_evoLandscape_gillespie_all.mat

white = 225;

Cmap = ([linspace(0,white,32) white*ones(1,1)  linspace(white,237,32);...
    linspace(114,white,32) white*ones(1,1)  linspace(white,177,32);...
    linspace(189,white,32) white*ones(1,1)  linspace(white,32,32)]'/255);

vec_xy = [1 -1; 1 0];

yv_vector = MODE.yv(MODE.yv>=-3 & MODE.yv<=-1);
xv_vector = MODE.xv(MODE.xv>=-1.5 & MODE.xv<=1.5);

%% plot evolutionary trajectories on fitness landscapes & fitness effects of promoter mutations

spx = 30;
spy = 80;

%which mutational probabilities?
m = 1;

for j=1:3
    
    if j==1
        plot_evolutionary_trajectories = [find(start_positions(:,1) == 41 & start_positions(:,2) == 2);
            find(start_positions(:,1) == 80 & start_positions(:,2) == 2);
            find(start_positions(:,1) == 80 & start_positions(:,2) == 30);
            find(start_positions(:,1) == 80 & start_positions(:,2) == 2);
            find(start_positions(:,1) == 61 & start_positions(:,2) == 2);]';
    elseif j==2
        plot_evolutionary_trajectories = [find(start_positions(:,1) == 41 & start_positions(:,2) == 60);
            find(start_positions(:,1) == 11 & start_positions(:,2) == 60);
            find(start_positions(:,1) == 80 & start_positions(:,2) == 30);
            find(start_positions(:,1) == 21 & start_positions(:,2) == 60);
            find(start_positions(:,1) == 61 & start_positions(:,2) == 60);]';
    elseif j==3
        plot_evolutionary_trajectories = [find(start_positions(:,1) == 41 & start_positions(:,2) == 2);
            find(start_positions(:,1) == 80 & start_positions(:,2) == 2);
            find(start_positions(:,1) == 80 & start_positions(:,2) == 30);
            find(start_positions(:,1) == 80 & start_positions(:,2) == 2);
            find(start_positions(:,1) == 41 & start_positions(:,2) == 60);]';
    end
    
    f=figure;
    
    F = principal_landscapes(:,:,1) * pc_combos(j,1) + principal_landscapes(:,:,2) * pc_combos(j,2);
    
    mut_x0 = 1-min([vec_xy(:,1); 0]):size(F,2)-max([vec_xy(:,1); 0]);
    mut_y0 = 1-min([vec_xy(:,2); 0]):size(F,1)-max([vec_xy(:,2); 0]);
    
    %%%%%%calculate fitness effects of mutations
    dF_freq = NaN(size(F));
    dF_freq(mut_y0+vec_xy(1,2),mut_x0+vec_xy(1,1)) = (F(mut_y0+vec_xy(1,2),mut_x0+vec_xy(1,1)) - F(mut_y0,mut_x0))/0.05; %whats the 0.05 for?
    dF_size = NaN(size(F));
    dF_size(mut_y0+vec_xy(2,2),mut_x0+vec_xy(2,1)) = (F(mut_y0+vec_xy(2,2),mut_x0+vec_xy(2,1)) - F(mut_y0,mut_x0))/0.05;
    
    %%%%%%plot fitness effects of burst size mutations
    s=subplot(3,2,1);
    Vy = linspace(-max(abs(quantile(abs(dF_size(:)),[.001 .999]))),max(abs(quantile(abs(dF_size(:)),[.001 .999]))),15);
    p=pcolor(landscape_windowed.X(mut_x0,1),yv_vector(mut_y0),dF_size(mut_y0,mut_x0));
    set(p, 'EdgeColor', 'none');
    hold on
    [cont,h]=contour(landscape_windowed.X(mut_x0,1),yv_vector(mut_y0),dF_size(mut_y0,mut_x0),Vy);
    %label iso fitness contour line
    set(h,'LineColor','black');
    h.LevelList = round(h.LevelList,3);
    clabel(cont,h,0,'LabelSpacing',25);
    caxis([-0.02 0.02])
    cbh=colorbar('v');
    set(cbh,'YTick',[-0.02 0 0.02])
    hold on
    axis([-1.5 1.5 -3 -1])
    set(s,'XTick',[-1 0 1],'YTick',[-3 -2 -1])
    title('\DeltaW(s^+)','FontSize',10)
    
    %find iso fitness contour line for later use
    idxc=1;
    idx=1;
    clear helper cont_value
    while idxc < length(cont)% && helper ~= 1
        helper{idx,1} = cont(:,idxc+1:idxc+cont(2,idxc));
        cont(1,idxc) = round(cont(1,idxc),3);
        cont_value(idx,1) = cont(1,idxc);
        idxc = idxc+cont(2,idxc)+1;
        idx = idx+1;
    end
    if ~isempty(helper(cont_value==0))
        dFy_cont10 = cell2mat(helper(cont_value==0)');
        dFy_cont10(:,~(dFy_cont10(1,:)>=-1 & dFy_cont10(1,:)<= 1.5) | ~(dFy_cont10(2,:)>=-3.1 & dFy_cont10(2,:)<= -1)) = [];
    else
        dFy_cont10 = [];
    end
    
    
    %%%%%%plot fitness effects of burst frequency mutations
    s=subplot(3,2,2);
    Vx = linspace(-max(abs(quantile(abs(dF_freq(:)),[.001 .999]))),max(abs(quantile(abs(dF_freq(:)),[.001 .999]))),15);
    p=pcolor(landscape_windowed.X(mut_x0,1),yv_vector(mut_y0),dF_freq(mut_y0,mut_x0));
    set(p, 'EdgeColor', 'none');
    hold on
    [cont,h]=contour(landscape_windowed.X(mut_x0,1),yv_vector(mut_y0),dF_freq(mut_y0,mut_x0),Vx);
    %label iso fitness contour line
    set(h,'LineColor','black');
    h.LevelList = round(h.LevelList,3);
    
    %     [cont,h]=contourf(landscape_windowed.X(mut_x0,1),yv_vector(mut_y0),dF_freq(mut_y0,mut_x0),Vx);
    %label iso fitness contour line
    clabel(cont,h,0,'LabelSpacing',25);
    caxis([-0.02 0.02])
    cbh=colorbar('v');
    set(cbh,'YTick',[-0.02 0 0.02])
    hold on
    axis([-1.5 1.5 -3 -1])
    set(s,'XTick',[-1 0 1],'YTick',[-3 -2 -1])
    title('\DeltaW(f^+)','FontSize',10)
    
    %find iso fitness contour line for later use
    idxc=1;
    idx=1;
    clear helper cont_value
    while idxc < length(cont)% && helper ~= 1
        helper{idx,1} = cont(:,idxc+1:idxc+cont(2,idxc));
        cont_value(idx,1) = round(cont(1,idxc),5);
        idxc = idxc+cont(2,idxc)+1;
        idx = idx+1;
    end
    if ~isempty(helper(cont_value==0))
        dFx_cont10 = cell2mat(helper(cont_value==0)');
        dFx_cont10(:,~(dFx_cont10(1,:)>=-1 & dFx_cont10(1,:)<= 1.5) | ~(dFx_cont10(2,:)>=-3.1 & dFx_cont10(2,:)<= -1)) = [];
    else
        dFx_cont10 = [];
    end
    
    
    %%%%%%plot evo. trajectory on fitness landscape
    cluster_col = linspecer(7);
    cluster_col_faint = cluster_col+(1-cluster_col)/3;
    cluster_assoc = [1 2 3 4];
    symbols = 'oxp*sdv^<+>h';
    
    
    s=subplot(3,3,[4 8]);
    
    %plot fitness landscape contour lines
    V = -.05:.005:.025;
    contour(landscape_windowed.X(:,1),yv_vector,F,V,'k')
    hold on
    %plot iso-fitness contour lines in size/frequency direction
    if ~isempty(dFx_cont10)
        plot(dFx_cont10(1,:),dFx_cont10(2,:),'k--')
    end
    if ~isempty(dFy_cont10)
        plot(dFy_cont10(1,:),dFy_cont10(2,:),'k-.')
    end
    
    %plot trajectories
    color_rand_idx = ceil(rand(size(evolutionary_path,2),1)*6);
    for i=plot_evolutionary_trajectories
        %trajectories
        p=plot(landscape_windowed.X(evolutionary_path{j,i,m}(:,2),1),yv_vector(evolutionary_path{j,i,m}(:,1)),...
            '-','Color',cluster_col(color_rand_idx(i),:),'MarkerSize',5,'LineWidth',.75);
        p.Color(4) = 0.5;
        %label starting points
        plot(landscape_windowed.X(evolutionary_path{j,i,m}(1,2),1),yv_vector(evolutionary_path{j,i,m}(1,1)),...
            'o','Color',cluster_col(color_rand_idx(i),:),'MarkerSize',5,'LineWidth',.75)
    end
    xlabel('expression rel. to WT [log2]')
    xlabel('mean expression [log2]')
    ylabel('noise [log2(CV)]')
    ylabel('noise [log2]')
    colormap(Cmap)
    caxis(V([1 end]))
    axis([-1.5 1.5 -3 -1])
    set(s,'XTick',[-1 0 1],'YTick',[-3 -2 -1])
    box off
    
    
    
    if j == 1
        print(f,'-depsc2','-painters','-loose',['results/plot_evotrajectories_PT1.eps'])
    elseif j==2
        print(f,'-depsc2','-painters','-loose',['results/plot_evotrajectories_PT2.eps'])
    elseif j==3
        print(f,'-depsc2','-painters','-loose',['results/plot_evotrajectories_peaked.eps'])
    end
    
    
    
    %%% time course comparision
    if j==3
        i=find(start_positions(:,1) == 80 & start_positions(:,2) == 2);
        i=i(1);
        
        f=figure;
        s=subplot(4,1,1);
        for m = [1 2 3]
            plot(evolutionary_path{j,i,m}(:,3),yv_vector(evolutionary_path{j,i,m}(:,1)),'Color',cluster_col(m,:))
            hold on
        end
        axis([0 2*10^6 -3 -1])
        set(s,'XTick',[0:0.5:2]*10^6,'XMinorTick','on','YTick',[-3 -2 -1],'FontSize',14)
        xlabel('time [a.u.]','FontSize',18)
        ylabel('noise [log2]','FontSize',18)
        box off
        legend('Pf=Ps','Pf=10*Ps','10*Pf=Ps')
        
        s=subplot(4,1,2);
        for m = [1 2 3]
            plot(evolutionary_path{j,i,m}(:,3),landscape_windowed.X(evolutionary_path{j,i,m}(:,2),1),'Color',cluster_col(m,:))
            hold on
        end
        axis([0 2*10^6 -1.5 1.5])
        set(s,'XTick',[0:0.5:2]*10^6,'XMinorTick','on','YTick',[-3 -2 -1],'FontSize',14)
        xlabel('time [a.u.]','FontSize',18)
        ylabel('noise [log2]','FontSize',18)
        box off
        
        s=subplot(4,3,[7 11]);
        %%%%%%plot evo. trajectory on fitness landscape
        cluster_col = linspecer(7);
        cluster_col_faint = cluster_col+(1-cluster_col)/3;
        cluster_assoc = [1 2 3 4];
        symbols = 'oxp*sdv^<+>h';
        
        F = principal_landscapes(:,:,1) * pc_combos(j,1) + principal_landscapes(:,:,2) * pc_combos(j,2);
        
        %plot fitness landscape contour lines
        V = -.05:.005:.025;
        contour(landscape_windowed.X(:,1),yv_vector,F,V,'k')
        hold on
        %plot iso-fitness contour lines in size/frequency direction
        if ~isempty(dFx_cont10)
            plot(dFx_cont10(1,:),dFx_cont10(2,:),'k--')
        end
        if ~isempty(dFy_cont10)
            plot(dFy_cont10(1,:),dFy_cont10(2,:),'k-.')
        end
        
        %plot trajectories
        color_rand_idx = 1:4;
        for m = [1 2 3]
            %trajectories
            p=plot(landscape_windowed.X(evolutionary_path{j,i,m}(:,2),1),yv_vector(evolutionary_path{j,i,m}(:,1)),...
                '-','Color',cluster_col(color_rand_idx(m),:),'MarkerSize',5,'LineWidth',.75);
            p.Color(4) = 0.5;
            %label starting points
            plot(landscape_windowed.X(evolutionary_path{j,i,m}(1,2),1),yv_vector(evolutionary_path{j,i,m}(1,1)),...
                'o','Color',cluster_col(color_rand_idx(m),:),'MarkerSize',5,'LineWidth',.75)
        end
        
        
        xlabel('mean expression [log2]','FontSize',18)
        
        ylabel('noise [log2]','FontSize',18)
        colormap(Cmap)
        caxis(V([1 end]))
        axis([-1.5 1.5 -3 -1])
        set(s,'XTick',[-1 0 1],'YTick',[-3 -2 -1],'FontSize',16)
        box off
        
        
        print(f,'-depsc2',['results/plot_trajectory_timecourse_peaked.eps'])
    end
    
end


%% epistasis between consecutive and opposing mutations

j = 3;
F = principal_landscapes(:,:,1) * pc_combos(j,1) + principal_landscapes(:,:,2) * pc_combos(j,2);

clear vec_xy2 F_ab
vec_xy2{1} = [1 -1; 1 -1];
vec_xy2{2} = [-1 0; -1 0];
vec_xy2{3} = [1 -1; -1 0];


%coordinates for double mutant epistasis plot
NM0 = [43 36];

f=figure;
%plot epistasis diagrams
for j=1
    s=subplot(2,2,4);
    
    for i=1:length(vec_xy2)
        
        
        mut_x2 = 2-min([vec_xy2{i}(:,1); 0]):(size(F,2)-2)-max([vec_xy2{i}(:,1); 0]);
        mut_y2 = 2-min([vec_xy2{i}(:,2); 0]):(size(F,1)-2)-max([vec_xy2{i}(:,2); 0]);
        
        
        dFx2 = (F(mut_y2+vec_xy2{i}(1,2),mut_x2+vec_xy2{i}(1,1)) - F(mut_y2,mut_x2));
        dFy2 = (F(mut_y2+vec_xy2{i}(2,2),mut_x2+vec_xy2{i}(2,1)) - F(mut_y2,mut_x2));
        
        dFxy2 = (F(mut_y2+sum(vec_xy2{i}(:,2)),mut_x2+sum(vec_xy2{i}(:,1))) - F(mut_y2,mut_x2));
        
        E = dFxy2 - dFx2 - dFy2;
        
        F_ab(1,1,i) = F(NM0(j,1),NM0(j,2));
        F_ab(1,2,i) = F(NM0(j,1)+vec_xy2{i}(2,2),NM0(j,2)+vec_xy2{i}(2,1));
        F_ab(2,1,i) = F(NM0(j,1)+vec_xy2{i}(1,2),NM0(j,2)+vec_xy2{i}(1,1));
        F_ab(2,2,i) = F(NM0(j,1)+sum(vec_xy2{i}(:,2)),NM0(j,2)+sum(vec_xy2{i}(:,1)));
        text(1+4*(i-1),F_ab(1,1,i),['g_0'],'Horizontal','right','Vertical','middle')
        hold on
        if i == 1
            text(2+4*(i-1),F_ab(2,1,i),['f^+'],'Horizontal','center','Vertical','middle')
            text(3+4*(i-1),F_ab(2,2,i),['f^+,f^+'],'Horizontal','left','Vertical','middle','Color','r')
        elseif i==2
            text(2+4*(i-1),F_ab(1,2,i),['s^-'],'Horizontal','center','Vertical','middle')
            text(3+4*(i-1),F_ab(2,2,i),['s^-,s^-'],'Horizontal','left','Vertical','middle','Color','r')
        elseif i==3
            text(2+4*(i-1),F_ab(1,2,i),['s^-'],'Horizontal','center','Vertical','middle')
            text(2+4*(i-1),F_ab(2,1,i),['f^+'],'Horizontal','center','Vertical','middle')
            text(3+4*(i-1),F_ab(2,2,i),['f^+,s^-'],'Horizontal','left','Vertical','middle','Color','r')
        end
        
        plot(4*(i-1)+[1.1 1.6],[F_ab(1,1,i) F_ab(1,2,i)],'k')
        
        plot(4*(i-1)+[1.1 1.6],[F_ab(1,1,i) F_ab(2,1,i)],'k')
        plot(4*(i-1)+[2.4 2.9],[F_ab(1,2,i) F_ab(2,2,i)],'r')
        plot(4*(i-1)+[2.4 2.9],[F_ab(2,1,i) F_ab(2,2,i)],'r')
        plot(4*(i-1)+[2.4 2.9],[F_ab(1,2,i) F_ab(2,1,i)+F_ab(1,2,i)-F_ab(1,1,i)],'k--')
        plot(4*(i-1)+[2.4 2.9],[F_ab(2,1,i) F_ab(2,1,i)+F_ab(1,2,i)-F_ab(1,1,i)],'k--')
    end
    F_ab_ex = [min(min(min(F_ab))) max(max(max(F_ab)))];
    
    
    s=subplot(2,2,3);
    
    [x,y]=meshgrid(landscape_windowed.X(mut_x2,1)+0.025,yv_vector(mut_y2));
    match = NaN(size(dFx2));
    match(dFx2<0 & dFy2>0 & E>0) = 3;
    match(dFx2>0 & dFy2>0 & E>0) = 1;
    match(dFx2>0 & dFy2<0 & E>0) = 2;
    match(dFx2<0 & dFy2<0 & E>0) = 4;
    
    plot(dFx_cont10(1,:),dFx_cont10(2,:),'k')
    hold on
    plot(dFy_cont10(1,:)-0.05,dFy_cont10(2,:),'k')
    
    contour(landscape_windowed.X(:,1),yv_vector,F,V,'k')
    text(landscape_windowed.X(NM0(1,2),1),yv_vector(NM0(1,1)),'1','Horizontal','center','Vertical','middle')
    
    title(['fitness, ' j],'FontSize',14)
    title('fitness','FontSize',12)
    
    xlabel('expression rel. to WT [log2]')
    xlabel('mean expression')
    ylabel('noise [log2(CV)]')
    ylabel('noise')
    colormap(Cmap)
    axis([-1.5 1.5 -3 -1])
    set(s,'XTick',[-1 0 1],'YTick',[-3 -2 -1])
    box off
    print(f,'-depsc2','-painters','-loose',['results/peaked_landscape_epistasis.eps'])
end