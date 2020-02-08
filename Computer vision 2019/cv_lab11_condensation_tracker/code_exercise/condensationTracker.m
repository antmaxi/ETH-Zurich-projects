function condensationTracker(videoName,params)
%condensationTracker(videoName,params)
%
% videoName  - videoName 
% params - parameters structure
%        . draw_plots {0,1} draw output plots throughout
%        . hist_bin   1-255 number of histogram bins for each color: proper values 4,8,16
%        . alpha      number in [0,1]; color histogram update parameter (0 = no update)
%        . sigma_position   std. dev. of system model position noise
%        . sigma_observe    std. dev. of observation model noise
%        . num_particles    number of particles
%        . model      {0,1} system model (0 = no motion, 1 = constant velocity)
%
% if using model = 1 then the following parameters are used:
%        . sigma_velocity   std. dev. of system model velocity noise
%        . initial_velocity initial velocity to set particles to
%
% Computer Vision - Autumn 2011
% Exercise 9 - Tracking
% Alex Mansfield and Bogdan Alexe

%%%oldAlpha = setFigTransparency(hFig,0.7);

% initialize the random generator
RandStream.setGlobalStream(RandStream('mt19937ar','seed',0));

% use AVI or WMV files?
use_wmv = true;

% load parameters
load ../data/params

% load video --------------------------------------------------------------
videoName = 'video3';
switch videoName
    case 'video1'
        % simple hand tracking on a pale background
        firstFrame = 10;
        lastFrame = 42;
        stepFrame = 1;
    case 'video2'
        % clutter and small occlusion
        firstFrame = 1;
        lastFrame = 40;
        stepFrame = 1;
    case 'video3'
        % non-constant velocity
        firstFrame = 1;
        lastFrame = 60;
        stepFrame = 1;
  case 'myOwnVideo'
    %implement here
end

% get the video
if(use_wmv)
	vid = VideoReader(['../data/' videoName '.wmv']); %  mmreader
else
    vid = aviread(['../data/' videoName '.avi']);
end

% get the first frame
if(use_wmv)
    frame = read(vid,firstFrame);
else
    frame = vid(firstFrame).cdata;
end
sizeFrame = size(frame);
heightFrame = sizeFrame(1);
widthFrame = sizeFrame(2);
frameValues = (firstFrame+stepFrame):stepFrame:lastFrame;

disp(sizeFrame);

% -------------------------------------------------------------------------

% tracking ----------------------------------------------------------------

figure(1);image(frame);

% USER INTERACTION
% draw initial bounding box and then double click
bb = imrect(gca);
wait(bb);
initialBB = round(getPosition(bb));

% bounding box size
WidthBB = initialBB(3);
HeightBB = initialBB(4);
disp(initialBB);
% GET INITIAL COLOR HISTOGRAM
%=== implement function color_histogram.m ===
%disp(initialBB(1));
%disp(params.hist_bin);
hist = color_histogram(initialBB(1),initialBB(2),initialBB(1)+initialBB(3),initialBB(2)+initialBB(4),frame,params.hist_bin);

%======================
params.model = 0;
disp(params.model);
state_length = 2;
if( params.model == 1 )
    state_length = 4;
end

% change velocity
params.initial_velocity = [5, 0];
params.sigma_velocity = 3;
params.sigma_position = 15;
params.sigma_observe = 0.05;
disp('velocity');
disp(params.initial_velocity);
meanStateAPriori = zeros(length(frameValues),state_length); % a priori mean state
meanStateAPosteriori = zeros(length(frameValues),state_length); % a posteriori mean state
meanStateAPriori(1,1:2) = [initialBB(1)+0.5*initialBB(3) initialBB(2)+0.5*initialBB(4)]; % bounding box centre

if (params.model==1 )    
    meanStateAPriori(1,3:4) = params.initial_velocity; % use initial velocity
end

% INITIALIZE PARTICLES
particles = repmat(meanStateAPriori(1,:), params.num_particles, 1);
particles_w = repmat(1/params.num_particles, params.num_particles, 1);
for i = 1:length(frameValues)
    disp(i)
    t = frameValues(i);
    
    % PROPAGATE PARTICLES
    %=== implement function propagate.m ===   
    tic
    disp('propagate')   
    particles = propagate(particles,sizeFrame,params);
    toc
    %======================
    
    % ESTIMATE
    %=== implement function estimate.m ===
    tic
    disp('estimate') 
    %disp(size(particles));
    meanStateAPriori(i,:) = estimate(particles, particles_w);
    disp('apriori')
    disp(meanStateAPriori(i,:));
    toc
    %======================
    
    % get frame
    if(use_wmv)
        frame = read(vid,t);
    else
        frame = vid(t).cdata;
    end
    
     % draw 
    if( params.draw_plots )
        figure(1);clf;image(frame);
        title(['Frame #' int2str(t)]);
        
        % plot a priori particles
        figure(1);hold on;        
        plot(particles(:,1),particles(:,2),'b.','MarkerSize',1);
        
        % plot a priori estimation
        figure(1);hold on;
        for j=i:-1:1
            lwidth = 30-3*(i-j);
            if(lwidth>0)                
                plot(meanStateAPriori(j,1),meanStateAPriori(j,2),'b.','MarkerSize',lwidth);
            end
            if(j~=i)                
                line([meanStateAPriori(j,1) meanStateAPriori(j+1,1)],[meanStateAPriori(j,2) meanStateAPriori(j+1,2)],'Color','b');
            end
        end
        
        % plot a priori bounding box
        if(~any(isnan(meanStateAPriori(i,:))))
            figure(1);hold on;           
            rectangle('Position',[meanStateAPriori(i,1)-0.5*WidthBB meanStateAPriori(i,2)-0.5*HeightBB WidthBB HeightBB],'EdgeColor','b');
        end
    end
    
    % OBSERVE
    %=== implement function observe.m ===
    tic
    disp('observe') 
    particles_w = observe(particles,frame,HeightBB,WidthBB,params.hist_bin,hist,params.sigma_observe);    
    toc
    %======================
    
    % UPDATE ESTIMATION  
    meanStateAPosteriori(i,:) = estimate(particles, particles_w);    
    disp('aposteriori')
    disp(meanStateAPosteriori(i,:));
    % update histogram color model
    tic
    disp('update');
    disp(min(max(1,round(meanStateAPosteriori(i,1)-0.5*WidthBB)),widthFrame));%, ...
    disp(min(max(1,round(meanStateAPosteriori(i,2)-0.5*HeightBB)),heightFrame));%, ...
    disp(min(max(1,round(meanStateAPosteriori(i,1)+0.5*WidthBB)),widthFrame));%, ...
    disp(min(max(1,round(meanStateAPosteriori(i,2)+0.5*HeightBB)),heightFrame));%
    hist_current = color_histogram(...
        min(max(1,round(meanStateAPosteriori(i,1)-0.5*WidthBB)),widthFrame), ...
        min(max(1,round(meanStateAPosteriori(i,2)-0.5*HeightBB)),heightFrame), ...
        min(max(1,round(meanStateAPosteriori(i,1)+0.5*WidthBB)),widthFrame), ...
        min(max(1,round(meanStateAPosteriori(i,2)+0.5*HeightBB)),heightFrame), ...
        frame,params.hist_bin);
    hist = (1-params.alpha).*hist + params.alpha.*hist_current;
    toc
    % draw 
    if( params.draw_plots )
        % plot weighted particles
        figure(1);hold on;        
        scatter(particles(:,1),particles(:,2),(eps+particles_w)*1e3,'r');
        
        % plot updated estimation
        figure(1);hold on;
        for j=i:-1:1
            lwidth = 30-3*(i-j);
            
            if(lwidth>0)                
                plot(meanStateAPosteriori(j,1),meanStateAPosteriori(j,2),'r.','MarkerSize',lwidth);
            end
            
            if(j~=i)                
                line([meanStateAPosteriori(j,1) meanStateAPosteriori(j+1,1)],[meanStateAPosteriori(j,2) meanStateAPosteriori(j+1,2)],'Color','r');
            end 
        end
        
        % plot updated bounding box
        if(~any(isnan(meanStateAPosteriori(i,:))))
            figure(1);hold on;            
            rectangle('Position',[meanStateAPosteriori(i,1)-0.5*WidthBB meanStateAPosteriori(i,2)-0.5*HeightBB WidthBB HeightBB],'EdgeColor','r');
        end
    end
    
    % RESAMPLE PARTICLES
    %=== implement function resample.m ===
    tic
    disp('resample');
    [particles, particles_w] = resample(particles,particles_w);
    toc
    %======================
      
    waitforbuttonpress;
    
end

% -----------------------------------------------------

end
