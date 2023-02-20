% paramPlot.m
% This script produces a parameter plot that shows all possible behaviors
% of the coupled Wilson-Cowan equations in WC4.
% For each point in a grid of values for gxy,gyx, it randomly produces a
% set of adjacency matrices, and for each one it obtains the corresponding
% 'behaviors' vector by calling ICsearch.

clear all; 
close all;

% ~~~~~~ parameter values ~~~~~~~
GENS = 5;  % maximum number of generations to iterate in ICsearch
c1 = 0.2; c2 = 0.3; % velocity update parameters in ICsearch
NICs = 20; % popluation size in ICsearch

tmin = 0;  % initial time of simulation
tmax = 80; % final time of simulation
N = 4;  % number of nodes per module
bx = 1.3;
by = 2;
thetax = 4;
thetay = 3.7;
gxx = 16/N; % = 4 -- strength of X self-connectivity
gyy = 3/N; % = 0.75 -- strength of Y self-connectivity
P = 1.5;  % external stimulus to excitatory units
%Q = 0; % This value isn't used

% creating gxy,gyx values
gN = 20;  % number of values to be evaluated for each of gxy,gyx
Gxy = linspace(0.01,30,gN);
Gyx = linspace(0.01,30,gN);

behavCount = zeros(gN*gN,6); % matrix specifying, for each (gxy,gyx) value, how
           % for how many adjacency matrices the behavior was observed.
           % The behaviors corresponding to the first four entries are the same
           % as the behaviors for the return value of 'ICsearch'. Namely, a
           % fixed point, multiple fixed points, periodic oscillations, and
           % non periodic oscillations.
           % The 5th entry corresponds to the combination of a
           % fixed point and oscillations. The sixth entry corresponds to
           % the combination of multiple fixed points and oscillations.
Nadj = 200;  % number of adjacency matrices to test

tic;
%~~~~~ MAIN LOOP ~~~~~
for xy_index = 1:gN
    for yx_index = 1:gN
        
        gxy = Gxy(xy_index);
        gyx = Gyx(yx_index);
        
        % sampling adjacency matrices
        for adj_index = 1:Nadj
            %densA = 8 + round(3*randn);  % number of 1 entries in A
            %densA = max(1,densA); densA = min(15,densA);
            %densB = 8 + round(3*randn);  % number of 1 entries in B
            %densB = max(1,densB); densB = min(15,densB);
            densA = 2;
            densB = 2;
            
            A = zeros(N*N,1);
            A(randperm(N*N,densA)) = 1;
            A = reshape(A,N,N);
            B = zeros(N*N,1);
            B(randperm(N*N,densB)) = 1;
            B = reshape(B,N,N);
            
            % some precomputed quantities
            Ayx = gyx*A;
            Bxy = gxy*B;
            fsx = 1/(1 + exp(bx*thetax));
            fsy = 1/(1 + exp(by*thetay));
            
            behav = ...
                ICsearch(N,bx,by,thetax,thetay,gxx,gyy,P,Ayx,Bxy,fsx,fsy,tmin,tmax,GENS,c1,c2,NICs);
            
            indy = gN*(xy_index-1)+yx_index;  % an index for behavCount
            
            behavCount(indy,1) = behavCount(indy,1) + (behav(1) > 0);
            
            behavCount(indy,2) = behavCount(indy,2) + (behav(2)-behav(1) > 0.07);
            
            behavCount(indy,3) = behavCount(indy,3) + (behav(3) > 0);
            
            behavCount(indy,4) = behavCount(indy,4) + (behav(4) > 0);
            
            behavCount(indy,5) = ...  
                behavCount(indy,5) + (behav(1)>0 && behav(3)>0);
            
            behavCount(indy,6) = ...  
                behavCount(indy,6) + ((behav(2)-behav(1)>0.07)&&(behav(3)>0));
     
        end
        disp(['Inner loop done with iteration ', num2str(yx_index)]);
    end
    disp(['>>> Outer loop done with iteration ', num2str(xy_index)]);
    %save('slow40c.mat','behavCount','gN','Gxy','Gyx','GENS','NICs','c1','c2','Nadj','tmin','tmax');
end
toc;

% plotting the results
for B = 1:6
    fig = zeros(gN,gN);   % the figure
    for idx1 = 1:gN
        for idx2 = 1:gN
            fig(idx1,idx2) = behavCount(gN*(idx1-1)+idx2,B);
        end
    end
    
    figure;
    imagesc(fig);
    colorbar;
end