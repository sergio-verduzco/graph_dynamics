function B = ICsearch(N,bx,by,thetax,thetay,gxx,gyy,P,Ayx,Bxy,fsx,fsy,tmin,tmax,GENS,c1,c2,NICs)
% B = ICsearch(...) receives the parameters of a WC4 system, along with the
% parameters for a particle swarm optimization search. For the given
% system, it does a search for behaviors of the WC4 system in the space of
% initial conditions.  Make sure NICs > 4.

% initial conditions
ICs = 0.1+zeros(2*N,NICs);
ICs(:,1) = repmat(0.05,2*N,1);
ICs(:,2) = [repmat(0.05,N,1);repmat(0.7,N,1)];
ICs(:,3) = [repmat(0.7,N,1);repmat(0.05,N,1)];
ICs(:,4) = repmat(0.7,2*N,1);
% the first 4 ICs are hand set, the rest are random
ICs(:,5:NICs) = rand(2*N,NICs-4);

% implementing swarm optimization
B = [0 0 0 0];   % each entry in this vector corresponds to a 
  % particular behavior. First: lowest fixed point. Second: highest fixed
  % point. Third: amplitude of periodic oscillation, if any. Fourth: max
  % amplitude of non-periodic oscillations.

F = zeros(NICs,3);  % fitness value for each initial condition
parfor i=1:NICs   % subsitute with 'for' to avoid parellel computing
    F(i,:) = WC4eval(ICs(:,i),N,bx,by,thetax,thetay,gxx,gyy,P,Ayx,Bxy,fsx,fsy,tmin,tmax);
end

if sum(F(:,1)==0) > 0  % if there are fixed points
    B(1) = min(F(F(:,1)==0, 2));   % smallest average value for a fixed point
    B(2) = max(F(F(:,1)==0, 2));   % largest average value for a fixed point
end
if sum(F(:,1)==1) > 0  % if there are periodic oscillations
    B(3) = max(F(F(:,1)==1, 3));   % largest oscillation amplitude
end
if sum(F(:,1)==2) > 0  % if there are non periodic oscillations
    B(4) = max(F(F(:,1)==2, 3));   % largest oscillation amplitude
end

%figure
%plot(F(:,3),'.');
util = 0; % util = 0; searching for fixed points
          % util = 1: searching for bistable fixed points 
          % util = 2: searching for periodic oscillations
          % util = 3: searching for non periodic oscillations
done = 0;  % done = 1 if all behaviors have been seen already

if B(1) == 0 && B(2) == 0
    util = 0;
else
    if B (2) - B(1) > 0.07  % bistability
        if B(3) > 0 % periodic oscillations too
            if B(4) > 0 % non periodic oscillations too
                done = 1;  % all behaviors found
            else
                util = 3;  % searching for non periodic oscillations
            end
        else
            util = 2;  % searching for periodic oscillations
        end
    else
        util = 1;
    end
end

vels = 0.05 - 0.1*rand(2*N,NICs);  % velocities. Same size as the ICs array

if util > 0   % if we're not searching for single fixed points
    [bestF, bestI] = max(F(:,3));  % max utility value
else
    [bestF, bestI] = min(F(:,3)); % best IC has the smallest oscillations
    bestF = 1 - bestF; % utility is 1 - (amplitud of smallest oscillation)
end
bestIC = ICs(:,bestI);   % initial condition with max utility value
curBestIC = bestIC;     % current IC with best utility value
%curBestF = bestF;   % current best utility

generation = 1;
while (generation  < GENS) && ~done
    %% UPDATE VELOCITIES AND ICs
    vels = vels + ...
          c1*rand*bsxfun(@minus,bestIC,ICs) + ...
          c2*rand*bsxfun(@minus,curBestIC,ICs);
    vels = min(vels,.4); % clipping fast velocities
    vels = max(vels,-0.4);
    
    ICs = ICs + vels;
    ICs = min(ICs,0.99);  % clipping large values
    ICs = max(ICs,0.001); % clipping small values
    
    %% EVALUATE ICs
    parfor i=1:NICs   % subsitute with 'for' to avoid parellel computing
        F(i,:) = WC4eval(ICs(:,i),N,bx,by,thetax,thetay,gxx,gyy,P,Ayx,Bxy,fsx,fsy,tmin,tmax);
    end
    
    %% UPDATE BEST ICs
    if util > 0   % if we're not searching for single fixed points
        [curBestF, curBestI] = max(F(:,3));  % max utility value
    else
        [curBestF, curBestI] = min(F(:,3)); % best IC has the smallest oscillations
        curBestF = 1 - curBestF;  % utility is 1 - (amplitud of smallest oscillation)
    end
    curBestIC = ICs(:,curBestI);   % IC with best utility from current generation
    if curBestF > bestF
        bestF = curBestF;
        bestIC = curBestIC;
    end
    
    %% UPDATE 'behaviors' VECTOR
    if sum(F(:,1)==0) > 0  % if there are fixed points
        B(1) = min(B(1), min(F(F(:,1)==0, 2))); 
        B(2) = max(B(2), max(F(F(:,1)==0, 2)));   
    end
    if sum(F(:,1)==1) > 0  % if there are periodic oscillations
        B(3) = max(B(3), max(F(F(:,1)==1, 3)));
    end
    if sum(F(:,1)==2) > 0  % if there are non periodic oscillations
        B(4) = max(B(4), max(F(F(:,1)==2, 3)));
    end
    
    %% SEE HOW MANY BEHAVIORS HAVE BEEN FOUND  
    if B(1) == 0 && B(2) == 0
        util = 0;
    else
        if B (2) - B(1) > 0.07  % bistability
            if B(3) > 0 % periodic oscillations too
                if B(4) > 0 % non periodic oscillations too
                    done = 1;  % all behaviors found
                else  % only missing non periodic oscillations
                    if util == 0
                        bestF = 0.001; % we are changing utility function, so bestF is reset
                    end
                    util = 3;  % searching for non periodic oscillations
                end
            else  % missing periodic oscillations
                if util == 0
                    bestF = 0.001; % we are changing utility function, so bestF is reset
                end
                util = 2;  % searching for periodic oscillations
            end
        else   % missing bistable fixed points
            if util == 0
                bestF = 0.001; % we are changing utility function, so bestF is reset
            end
            util = 1;
        end
    end
    
    generation = generation + 1;
end

