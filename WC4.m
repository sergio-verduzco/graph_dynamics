function Y = WC4(t,X,N,bx,by,thetax,thetay,gxx,gyy,P,Ayx,Bxy,fsx,fsy)
% Y = WC4(t,X) implements the right-hand side of the dynamic equations for
% the Wilson-Cowan system with two fully-connected modules, 4 nodes per
% module.
% The argument t is an unused time variable.
% The argument X is a vector with 8 elements denoting the state of the
% system. The first 4 elements correspond to the 'x' variables in Eq. 1 of
% the paper, whereas the next 4 elements correspond to the 'y' variables.
% The output Y corresponds to the time derivatives of the 8 elements in X.

% The matrices Ayx and Bxy are not the A and B binary matrices from
% equation 1. Ayx = A*g_{yx}, and B = B*g_{xy}. This avoids the extra
% multiplication.
% fsx = 1/(1 + exp(bx*thetax)). Included as a variable to save operations.
% Notice there's no Q argument, since we always use Q=0.


Xk = X(1:N); % the four x_k values of Eq. 1
Yk = X(N+1:2*N); % the four x_y values of Eq. 1

argX = -Ayx*Yk + P + gxx*sum(Xk);  % arguments to the x_k sigmoidals
argY =  Bxy*Xk + gyy*sum(Yk);      % arguments to the y_k sigmoidals

SX = ones(N,1)./(1 + exp(-bx*(argX - thetax))) - fsx;
SY = ones(N,1)./(1 + exp(-by*(argY - thetay))) - fsy;

Y = zeros(2*N,1);
Y(1:N) = -Xk + (1-Xk).*SX;
Y(N+1:2*N) = -Yk + (1-Yk).*SY;