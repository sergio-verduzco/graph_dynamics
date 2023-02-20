% densPlot.m
% This script produces a plot that shows all possible behaviors
% of the coupled Wilson-Cowan equations in WC4 as the density in the 
% adjacency matrices is changed.

clear all; 
close all;

% ~~~~~~ parameter values ~~~~~~~
GENS = 2;  % maximum number of generations to iterate in ICsearch
c1 = 0.2; c2 = 0.3; % velocity update parameters in ICsearch
NICs = 5; % popluation size in ICsearch

tmin = 0;  % initial time of simulation
tmax = 60; % final time of simulation
N = 4;  % number of nodes per module
bx = 1.3;
by = 2;
thetax = 4;
thetay = 3.7;
gxx = 16/N; % = 4 -- strength of X self-connectivity
gyy = 3/N; % = 0.75 -- strength of Y self-connectivity
gxy = 15;
gyx = 15;
P = 1.5;  % external stimulus to excitatory units
%Q = 0; % This value isn't used

densesA =6:10;  % vector of densities to be tested for matrix A
densesB = 6:8;  % vector of densities to be tested for matrix B
nA = length(densesA); nB = length(densesB);

behavCount = zeros(nA*nB,6); 
           % matrix specifying, for each density combination, how
           % many adjacency matrices with the behavior were observed.
           % The behaviors corresponding to the first four entries are the same
           % as the behaviors for the return value of 'ICsearch'. Namely, a
           % fixed point, multiple fixed points, periodic oscillations, and
           % non periodic oscillations.
           % The 5th entry corresponds to the combination of a
           % fixed point and oscillations. The sixth entry corresponds to
           % the combination of multiple fixed points and oscillations.
Nadj = 20;  % number of adjacency matrices to test for each density

tic;
%~~~~~ MAIN LOOP ~~~~~
for densA = densesA
    for densB = densesB
        % sampling adjacency matrices
        maxAdjA = nchoosek(N*N,densA);
        maxAdjB = nchoosek(N*N,densB);
        
        for adj_index = 1:min([Nadj maxAdjA maxAdjB]);
            
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
            
            indy = nB*(densA-densesA(1))+densB-densesB(1)+1;  % an index for behavCount
            
            behavCount(indy,1) = behavCount(indy,1) + (behav(1) > 0);
            
            behavCount(indy,2) = behavCount(indy,2) + (behav(2)-behav(1) > 0.07);
            
            behavCount(indy,3) = behavCount(indy,3) + (behav(3) > 0);
            
            behavCount(indy,4) = behavCount(indy,4) + (behav(4) > 0);
            
            behavCount(indy,5) = ...  
                behavCount(indy,5) + (behav(1)>0 && behav(3)>0);
            
            behavCount(indy,6) = ...  
                behavCount(indy,6) + ((behav(2)-behav(1)>0.07)&&(behav(3)>0));
     
        end
        disp(['Inner loop done with iteration ', num2str(densB)]);
    end
    disp(['>>> Outer loop done with iteration ', num2str(densA)]);
    %save('xy15yx15_densA5-12.mat','behavCount','N','densesA','densesB','GENS','NICs','c1','c2','Nadj','tmin','tmax');
end
toc;

% plotting the results
for B = 1:6
    fig = zeros(nA,nB);   % the figure
    for idx1 = 1:nA
        for idx2 = 1:nB
            fig(idx1,idx2) = behavCount(nB*(idx1-1)+idx2,B);
        end
    end
    
    figure;
    imagesc(fig);
    colorbar;
end