function F = WC4eval(X0,N,bx,by,thetax,thetay,gxx,gyy,P,Ayx,Bxy,fsx,fsy,tmin,tmax)
% F = WC4fitness1(X0,tmin,tmax)
% Given initial conditions X0, this function returns the type of dynamic
% behavior in F(1), the average activity for the first node in F(2), and a 
% measure of fitness in F(3) for the WC4 system that starts at X0 from time
% tmin to time tmax. The current fitness measure is the maximum difference
% between the mean value of a node activity and it's maximum deviation from
% that mean value.
% F(1) = 0 --> fixed point
% F(1) = 1 --> periodic oscillations
% F(1) = 2 --> non periodic oscillations

% To avoid false detection of non periodic oscillations arising from
% transients, whenever a isperiodic returns 0 the simulation will be run
% for a longer time to allow possible transients to fade.

[T,Y] = ode45(@(t,X)WC4(t,X,N,bx,by,thetax,thetay,gxx,gyy,P,Ayx,Bxy,fsx,fsy),[tmin tmax],X0);

tmid = tmin + (tmax-tmin)/2;
resY = Y(T>tmid,:);  
resT = T(T>tmid);
avg = mean(resY);
mx = max(resY(resT>(tmin + 0.75*(tmax-tmin)),:));
width = mx - avg;
[max_width, max_ind] = max(width);

F = zeros(3,1);
F(2) = avg(1);
F(3) = max_width;
if F(3) > 0.05   % oscillations
    F(1) = 2 - isperiodic(resT,resY(:,max_ind));  
else
    F(1) = 0;  % fixed point
end

if F(1) == 2  % if non periodic oscillations were detected
    X0 = Y(end,:)';
    tmax = 2*tmax;
    [T,Y] = ode45(@(t,X)WC4(t,X,N,bx,by,thetax,thetay,gxx,gyy,P,Ayx,Bxy,fsx,fsy),[tmin tmax],X0);
    
    tmid = tmin + (tmax-tmin)/2;
    resY = Y(T>tmid,:);
    resT = T(T>tmid);
    avg = mean(resY);
    mx = max(resY(resT>(tmin + 0.75*(tmax-tmin)),:));
    width = mx - avg;
    [max_width, max_ind] = max(width);
    
    F(2) = avg(1);
    F(3) = max_width;
    if F(3) > 0.05   % oscillations
        F(1) = 2 - isperiodic(resT,resY(:,max_ind));
    else
        F(1) = 0;  % fixed point
    end
end