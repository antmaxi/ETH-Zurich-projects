function matchingCost = shape_matching(X,Y,display_flag)

%computes the matching cost between template contour points X and target contour points Y  

%%%
%%%Define flags and parameters:
%%%

if nargin < 3
    display_flag = 0;
end

nbSamples = 100;%size(X,1);
nbBins_theta = 12;
nbBins_r = 5;
smallest_r = 1/8;%length of the smallest radius (assuming normalized distances)
biggest_r = 3;%length of the biggest radius (assuming normalized distances)
maxIterations = 3;

n_out = 0; % number of "dummy" points
lambda_K = 1;% how much times lambda is bigger than mean_distance^2
%%%%%%%%%%%%%  MODE FLAGS
tangent_flag = 0; % if true => to use tangent for angle descriptors
local_flag = 0; % if true => to use angle between direction to point and X axis,
                    % else angle of rotation around center of image
if tangent_flag
    local_flag = 1;
end
%pause on % to enable pause function
if display_flag
  % for i = 1:size(X,1)
       subplot(1,2,1)
       plot(X(:,1),X(:,2),'b+');%, X(i,1),X(i,2),'ro')
       axis('ij'), title('X');
       subplot(1,2,2)
       plot(Y(:,1),Y(:,2),'ro');%, Y(i,1),Y(i,2),'b+')
       axis('ij'), title('Y');
       drawnow	
       %pause(0.05);
  % end
end

if display_flag
   [x,y] = meshgrid(linspace(min(X(:,1))-10,max(X(:,1))+10,36),...
                    linspace(min(X(:,2))-10,max(X(:,2))+10,36));
   x = x(:);
   y = y(:);
   M = length(x);
end

%%%
%%% compute correspondences
%%%
X = datasample(X, nbSamples, 'Replace', false);

Y = datasample(Y, nbSamples, 'Replace', false);
% currentY = [Y; Inf*ones(n_out,2)];
% currentX = [X; Inf*ones(n_out,2)];
% currentX_pure = X; %without dummy points
currentIteration = 1;
currentX = X;
disp('computing shape contexts for target...')
ShapeDescriptors2 = sc_compute(Y',nbBins_theta,nbBins_r, ...
        smallest_r,biggest_r, tangent_flag, local_flag, n_out);
disp('done.')

deleted = 0;
while currentIteration <= maxIterations 
   % add dummy points to X and Y
   if n_out
        %X = [X; Inf*ones(n_out,2)];
        if currentIteration ~= 1
            %X = [X; %currentX = [X; datasample(X, n_out, 'Replace', false)];
        end
   end
   disp(['iter=' int2str(currentIteration)]);

   %write the sc_compute.m function
   disp('computing shape contexts for (deformed) model...')
   ShapeDescriptors1 = sc_compute(currentX',nbBins_theta,nbBins_r, ...
        smallest_r,biggest_r, tangent_flag, local_flag, deleted);
   disp('done.')

   %set lambda here
   mean_dist = mean(sqrt(dist2(currentX,currentX)), 'all');
   %disp('mean dist');
   %disp(mean_dist);
   lambda = lambda_K*mean_dist^2;  
   %write the chi2_cost.m function
   costMatrixC = chi2_cost(ShapeDescriptors1, ShapeDescriptors2, n_out, tangent_flag);
   disp(size(ShapeDescriptors1));
   disp(size(ShapeDescriptors2));
   
   % assigning labels
   corespondencesIndex = hungarian(costMatrixC);
   
%    nx = size(currentX, 1);
%    ny = size(Y, 1);
%    index = corespondencesIndex(corespondencesIndex <= ny-n_out);
%    mask_out = corespondencesIndex(corespondencesIndex >  ny - n_out);
%    mask_in = corespondencesIndex(corespondencesIndex <=  ny - n_out);
    Xwarped = currentX(corespondencesIndex,:);   
    Xunwarped = X(corespondencesIndex,:);
%    
%    % deleting outliers
%    Xunwarped = Xunwarped(1:nx - n_out, :);
%    Xwarped = Xwarped(1:nx - n_out, :);
%    X = X(1:nx - n_out, :);  
%    Y = Y(1:ny - n_out,:);
   
   if display_flag       
      figure(2)
      plot(Xwarped(:,1),Xwarped(:,2),'b+',Y(:,1),Y(:,2),'ro');
      hold on
      plot([Xwarped(:,1) Y(:,1)]',[Xwarped(:,2) Y(:,2)]','k-');      
      hold off
      axis('ij')
      title('correspondences (warped X)')
      drawnow	
   end
   
   %Xwarped(mask_out) = [];
   %disp(mask_out);
   %Xunwarped(mask_out) = [];
   
   if display_flag
      % show the correspondences between the untransformed images
      figure(3)
      plot(X(:,1),X(:,2),'b+',Y(:,1),Y(:,2),'ro')      
      hold on
      plot([Xunwarped(:,1) Y(:,1)]',[Xunwarped(:,2) Y(:,2)]','k-')
      hold off
      axis('ij')
      title(' correspondences (unwarped X)')
      drawnow	
   end
  
   [w_x,w_y,E] = tps_model(Xunwarped,Y,lambda);
   
   % warp each coordinate
   eps = 1e-10;
   fx_aff = w_x(nbSamples+1:nbSamples+3)'*[ones(1,nbSamples); X'];
   d2 = max(dist2(Xunwarped,X),0);  
   U = d2.*log(d2+eps);
   %disp(U);
   fx_wrp = w_x(1:nbSamples)'*U;
   fx = fx_aff+fx_wrp;   
   
   fy_aff = w_y(nbSamples+1:nbSamples+3)'*[ones(1,nbSamples); X'];
   fy_wrp = w_y(1:nbSamples)'*U;
   fy = fy_aff+fy_wrp;

   % update currentX for the next iteration
   currentX = [fx; fy]';   

   if display_flag
      figure(4)
      plot(currentX(:,1),currentX(:,2),'b+',Y(:,1),Y(:,2),'ro');
      axis('ij')
      title(['iteration=' int2str(currentIteration) ',  I_f=' num2str(E) ])
      % show warped coordinate grid
      fx_aff=w_x(nbSamples+1:nbSamples+3)'*[ones(1,M); x'; y'];
      d2 = dist2(Xunwarped,[x y]);
      fx_wrp = w_x(1:nbSamples)'*(d2.*log(d2+eps));
      fx = fx_aff+fx_wrp;
      fy_aff = w_y(nbSamples+1:nbSamples+3)'*[ones(1,M); x'; y'];
      fy_wrp = w_y(1:nbSamples)'*(d2.*log(d2+eps));
      fy = fy_aff+fy_wrp;
      hold on
      plot(fx,fy,'k.','markersize',1)
      hold off
      drawnow
   end         

   %update currentIteration
   currentIteration = currentIteration + 1;
   
end

matchingCost = E;