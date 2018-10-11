function fitness_sensitivity_noise_mean_changes()
% analyse sensitivity of fitness to changes in noise or mean expression on cluster-average fitness landscapes
%%% df/dx and df/dy function of fitness contours %%%

load data/FNM_smooth_fun_grid05_025.mat
load data/FNM_smooth_fun_grid05_025_landscape_PCA.mat


V = -.05:.005:.025;

Col2 = (linspecer(numel(V)));
%markersize
ms = 4;

pc_combos = [1 0;
    0 1;
    1 1];

pc_combo_labels = {'PC 1','PC 2','PC1 + PC2'};

%%

for u=1:size(pc_combos,1)
    
    
    PClandscape = pc_combos(u,1) * principal_landscapes(:,:,1) + pc_combos(u,2) * principal_landscapes(:,:,2);
    [PClandscape_dx,PClandscape_dy] = gradient(PClandscape,diff(MODE.xv([1 2])),diff(MODE.yv([1 2])));
    
    %decode contour object
    idxc=1;
    idx=1;
    idxl=1;
    clear cont_cell cont_values v L legendlabel cont_indicies
    cont = contourc(landscape_windowed.X(:,1),MODE.yv(MODE.yv >= -3 & MODE.yv <= -1),...
        PClandscape,V);
    
    while idxc < length(cont)
        
        h = cont(:,idxc+1:idxc+cont(2,idxc));
        %gating
        h(:,h(1,:) < -1.5 | h(1,:) > 1.5) = [];
        h(:,h(2,:) < -3 | h(2,:) > -1) = [];
        [~,helper] = min(abs(V-cont(1,idxc)));
        if ~isempty(h) && V(helper) >= -0.033
            [~,v(idx,1)] = min(abs(V-cont(1,idxc)));
            %approximate expr rel to wildtype
            [Xv,Xh]=meshgrid(landscape_windowed.X(:,1),h(1,:));
            [~,mindistX] = min(abs(Xv-Xh),[],2);
            %approximate noise
            [Yv,Yh]=meshgrid(MODE.yv(MODE.yv >= -3 & MODE.yv <= -1),h(2,:));
            [~,mindistY] = min(abs(Yv-Yh),[],2);
            cont_indicies{idx} = [mindistX mindistY];
            cont_values{idx} = h';
            idx=idx+1;
        end
        idxc = idxc+cont(2,idxc)+1;
    end
    
    if mod(u,3)==1
        f=figure;
    end
    
    
    %1) plot contour lines
    s=subplot(3,3,u);
    [c,h] = contour(landscape_windowed.X(:,1),MODE.yv(MODE.yv >= -3 & MODE.yv <= -1),...
        PClandscape,V);
    colormap(Col2)
    caxis([V(1) V(end)])
    
    hold on
    for i=1:numel(v)
        [~,max_noise] = max(cont_values{i}(:,2));
        plot(cont_values{i}(max_noise,1),cont_values{i}(max_noise,2),...
            's','Color',Col2(v(i),:),'MarkerFaceColor',Col2(v(i),:),'MarkerSize',ms)
        [~,max_mean] = max(cont_values{i}(:,1));
        plot(cont_values{i}(max_mean,1),cont_values{i}(max_mean,2),...
            'x','Color',Col2(v(i),:),'MarkerFaceColor',Col2(v(i),:),'MarkerSize',ms)
        [~,min_mean] = min(cont_values{i}(:,1));
        plot(cont_values{i}(min_mean,1),cont_values{i}(min_mean,2),...
            'v','Color',Col2(v(i),:),'MarkerFaceColor',Col2(v(i),:),'MarkerSize',ms)
        nr = 15;
        idx = 3:nr:numel(cont_values{i})/2;
        if numel(idx)>0
            dx = diag(nanmean(PClandscape_dx(cont_indicies{i}(idx,2),cont_indicies{i}(idx,1)),3))*10;
            dy = diag(nanmean(PClandscape_dy(cont_indicies{i}(idx,2),cont_indicies{i}(idx,1)),3))*10;
            quiver(cont_values{i}(idx,1)-dx/2,cont_values{i}(idx,2)-dy/2,...
                dx,dy,...
                0,'Color',[0.5 0.5 0.5])
        end
        
    end
    axis([-1.5 1.5 -3 -1])
    
    title(pc_combo_labels(u),'FontSize',8)
    ylabel('noise','FontSize',8)
    xlabel('mean expression','FontSize',8)
    box off
    set(s,'XTick',[-1:1],'YTick',[-3:-1],'FontSize',8)
    
    
    
    %2) plot df/dx and df/dy along contour lines
    s=subplot(3,3,u+3);
    limX = [-0.025 0.025];
    limY = [-0.029 0];
    for i=1:numel(v)
        plot(diag(nanmean(PClandscape_dx(cont_indicies{i}(:,2),cont_indicies{i}(:,1)),3)),...
            diag(nanmean(PClandscape_dy(cont_indicies{i}(:,2),cont_indicies{i}(:,1)),3)),...
            '-','Color',Col2(v(i),:),'LineWidth',1);
        
        limX(1) = min([limX(1); diag(nanmean(PClandscape_dx(cont_indicies{i}(:,2),cont_indicies{i}(:,1)),3))]);
        limX(2) = max([limX(2); diag(nanmean(PClandscape_dx(cont_indicies{i}(:,2),cont_indicies{i}(:,1)),3))]);
        limY(1) = min([limY(1); diag(nanmean(PClandscape_dy(cont_indicies{i}(:,2),cont_indicies{i}(:,1)),3))]);
        limY(2) = max([limY(2); diag(nanmean(PClandscape_dy(cont_indicies{i}(:,2),cont_indicies{i}(:,1)),3))]);
        
        hold on
        [h,max_noise] = max(cont_values{i}(:,2));
        plot(diag(nanmean(PClandscape_dx(cont_indicies{i}(max_noise,2),cont_indicies{i}(max_noise,1)),3)),...
            diag(nanmean(PClandscape_dy(cont_indicies{i}(max_noise,2),cont_indicies{i}(max_noise,1)),3)),...
            's','Color',Col2(v(i),:),'MarkerFaceColor',Col2(v(i),:),'MarkerSize',ms);
        
        [~,max_mean] = max(cont_indicies{i}(:,1));
        plot(diag(nanmean(PClandscape_dx(cont_indicies{i}(max_mean,2),cont_indicies{i}(max_mean,1)),3)),...
            diag(nanmean(PClandscape_dy(cont_indicies{i}(max_mean,2),cont_indicies{i}(max_mean,1)),3)),...
            'x','Color',Col2(v(i),:),'MarkerFaceColor',Col2(v(i),:),'MarkerSize',ms);
        [~,min_mean] = min(cont_indicies{i}(:,1));
        plot(diag(nanmean(PClandscape_dx(cont_indicies{i}(min_mean,2),cont_indicies{i}(min_mean,1)),3)),...
            diag(nanmean(PClandscape_dy(cont_indicies{i}(min_mean,2),cont_indicies{i}(min_mean,1)),3)),...
            'v','Color',Col2(v(i),:),'MarkerFaceColor',Col2(v(i),:),'MarkerSize',ms);
    end
    
    %plot reference lines
    plot(zeros(2,1),[limY(1)-0.001 limY(2)+0.001],'k--')
    plot([limX(1)-0.005 limX(2)+0.001],zeros(2,1),'k--')
    r=refline(-1,0);
    set(r,'Color','k','LineStyle',':')
    r=refline(1,0);
    set(r,'Color','k','LineStyle',':')
    
    
    
    axis([limX(1)-0.001 limX(2)+0.001 limY(1)-0.0005 limY(2)+0.0005])
    xlabel('fitness sensitivity to mean changes (x100)','FontSize',8)
    ylabel('fitness sensitivity to noise changes (x100)','FontSize',8)
    
    box off
    set(s,'XTick',[-.03:.01:0.03],'XTickLabel',[-.03:.01:0.03]*100,'YTick',[-.03:.01:0.03],'YTickLabel',[-.03:.01:0.03]*100,'FontSize',8)
end
print(f,'-depsc2','-painters','-loose',['results/fitness_sensitivities_PT_dfdx_dfdy.eps'])


%% sensitivities across all landscapes
clear lw_dx lw_dy


for u=1:size(pc_combos,1)
    PClandscape = pc_combos(u,1) * principal_landscapes(:,:,1) + pc_combos(u,2) * principal_landscapes(:,:,2);
    [PClandscape_dx(:,:,u),PClandscape_dy(:,:,u)] = gradient(PClandscape,diff(MODE.xv([1 2])),diff(MODE.yv([1 2])));
    
    PCmean_dx(u) = mean(mean(abs(PClandscape_dx(:,:,u))));
    PCmean_dy(u) = mean(mean(abs(PClandscape_dy(:,:,u))));
    
    PCmean_dx_m_dy(u) = mean(mean(abs(PClandscape_dx(:,:,u)) - abs(PClandscape_dy(:,:,u))));
    PCmean_dx_p_dy(u) = mean(mean(abs(PClandscape_dx(:,:,u)) + abs(PClandscape_dy(:,:,u))));
end

for i=1:33
    [lw_dx(:,:,i),lw_dy(:,:,i)] = gradient(landscape_windowed.F2d(MODE.yv>=-3 & MODE.yv<=-1,:,i),diff(MODE.xv([1 2])),diff(MODE.yv([1 2])));
    
    mean_dx(i) = mean(mean(abs(lw_dx(:,:,i))));
    mean_dy(i) = mean(mean(abs(lw_dy(:,:,i))));
    
    mean_dx_m_dy(i) = mean(mean(abs(lw_dx(:,:,i)) - abs(lw_dy(:,:,i))));
    mean_dx_p_dy(i) = mean(mean(abs(lw_dx(:,:,i)) + abs(lw_dy(:,:,i))));
end

f=figure;
scatter(mean_dx_p_dy,mean_dx_m_dy,50,'k','filled'),
hold on,
scatter(PCmean_dx_p_dy,PCmean_dx_m_dy,100,'r','filled')
plot([0,0.06],[0 0],'k--')
plot([0,0.04],[0 0.04],'k:')
plot([0,0.04],[0 -0.04],'k:')
axis([0 0.06 -0.04 0.04])
xlabel('mean combined sensitivities')
ylabel('mean differential sensitivities (mean - noise)')
legend('individual landscapes','principal topologies')
for u=1:size(pc_combos,1)
    text(PCmean_dx_p_dy(u),PCmean_dx_m_dy(u),num2str(u))
end
print(f,'-depsc2','-painters','-loose',['results/fitness_sensitivities_andscapes_dfdx_dfdy.eps'])
