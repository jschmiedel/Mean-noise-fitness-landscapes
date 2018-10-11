function Z=fitness_noise_mean_smoothing_v5(x,y,z,err_x,err_y,err_z,MODE)

%calculate z (fitness) values on grid points based on multi variate normal
%distributions taking errors in x,y and z into account

%either supply gridXY, a set of (irregular) points where to evaluate the fitness
%landscape, or the marginals of a grid "xv" and "yv" to let the script
%calculate a grid

% the smaller dist_weight, the smoother the surface
if ~isfield(MODE,'weight_xy')
    MODE.weight_xy = iqr(x)/iqr(y);
end

if ~isfield(MODE,'gridXY')
    [X,Y] = meshgrid(MODE.xv,MODE.yv);
    gridXY = [X(:) Y(:)];
else
    gridXY = MODE.gridXY;
end

if ~isfield(MODE,'dist_weight_xy')
    MODE.dist_weight_xy = [1 1];
elseif numel(MODE.dist_weight_xy) == 1
    MODE.dist_weight_xy = ones(1,2)*MODE.dist_weight_xy;
end 

%calculate weights based on mvn distributions
% 1) define nec. vectors
%create local surrounding of gridpoint
[Xlocal,Ylocal] = meshgrid(linspace(-3,3,MODE.local_grid_n),linspace(-3,3,MODE.local_grid_n));
localXY = [Xlocal(:) Ylocal(:)];

localXY(:,1) = localXY(:,1)*sqrt(MODE.dist_weight_xy(1));
localXY(:,2) = localXY(:,2)*sqrt(MODE.dist_weight_xy(2));
localXY_mvnpdf = mvnpdf(localXY,0,MODE.dist_weight_xy);

localXY = localXY(localXY_mvnpdf/max(localXY_mvnpdf) > 0.01,:);
localXY_mvnpdf = localXY_mvnpdf(localXY_mvnpdf/max(localXY_mvnpdf) > 0.01);

%promoter positions
xy = [x y];
%and error
err_xy = [err_x err_y];

% tic
% 2) weights of each promoter for a single grid point
W = cell2mat(arrayfun(@(i) (localXY_mvnpdf'*...
    reshape(mvnpdf(repmat(repmat(gridXY(i,:),size(localXY,1),1) + localXY,size(xy,1),1),...
    kron(xy,ones(size(localXY,1),1)),...
    permute(kron(err_xy,ones(size(localXY,1),1)),[3 2 1])),...
    size(localXY,1),size(xy,1))./...
    (err_z.^2)')',...
    1:size(gridXY,1),'uni',false));
% toc

%smoothed fitness at all grid points
Z = W'*z./sum(W,1)';

%if grid was calculated in script, reshape smoothed fitness to adhere to
%grid layout and remove grid points to far from actual data
if ~isfield(MODE,'gridXY')
    Z = reshape(Z,numel(MODE.yv),numel(MODE.xv));
    %discard grid points too far away from data
    if MODE.Xtra_range~=inf
        d = pdist2([X(:) Y(:)*MODE.weight_xy],[x y*MODE.weight_xy],'euclidean');
        Dmin = reshape(min(d,[],2),numel(MODE.yv),numel(MODE.xv));
        Z(Dmin>MODE.Xtra_range) = NaN;
    end
end

