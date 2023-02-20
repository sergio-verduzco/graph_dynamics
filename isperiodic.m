function P = isperiodic(T,Y)
% P = isperiodic(T,Y) receives a function of time specified by the vectors
% T,Y, where T(i) is the time when the function adopts the value Y(i). If
% the function is (roughly) periodic, then P=1; otherwise P=0.

% 1) Discretize the domain of function
N = 200; % number of points in domain
times = linspace(T(1),T(end),N);
Ydisc = interp1(T,Y,times);
Ydisc = Ydisc - mean(Ydisc);  % removing 'DC' component

% If you want to see the discretized, normalized signal
% figure;
% plot(T,Y);
% hold on
% plot(times,Ydisc,'r');

% 2) Calculate the normalized "convolution" of the function with itself
%   (it's actually the convolution of Y(i) with Y(-i)) and see if it is
%   larger than the set threshold
P = 0;  % default return value
thresh = 0.95;  % periodic if the convolution surpasses this value
pad = ceil(N/5); % number of points in the first sum
convY = zeros(1,N-2*pad+1);  % this vector will have the convolution
integY1 = sum(Ydisc(1:pad-1).*Ydisc(1:pad-1)); % integral of Y squared up to the point where the 
%                                              % convolution is being calculated                                             
integY2 = sum(Ydisc(N-pad+2:N).*Ydisc(N-pad+2:N)); % same, but from the back

for i = pad:N-pad;
    integY1 = integY1 + Ydisc(i)*Ydisc(i);
    integY2 = integY2 + Ydisc(N-i+1)*Ydisc(N-i+1);
    convY(i-pad+1) = (sum(Ydisc(1:i).*Ydisc(N-i+1:N)))/max(integY1,integY2);
    
    if convY(i-pad+1) > thresh
        P = 1;
        break;
    end
end

% figure
% plot(convY);

