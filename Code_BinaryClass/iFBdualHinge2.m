function [u,U,ww,W,wn,M,time]=iFBdualHinge2(w0,X,alpha,gamma,lambda0,T)

d=length(w0); n=length(X(:,1)); 

%Initialization
w = w0 ; u=X*w0 ; z=u;  U=u.*ones(n,T); W=w0.*ones(d,T);   % w0 primal, u0 dual, dual matrix, primal matrix
wn=w0.*ones(d,T);          % normalized error from solution wnast
A=X*(X');             %Gradient
%D = ( (norm(-X'*u)^2)/2 + sum(lambda0*u)/lambda0 ).*ones(T,1) ;   %Dual objective
M=min(X*w0).*ones(T,1);          %Margin                 
time=zeros(T-1,1) ;

for t=2:T
    tic ;
    lambda=lambda0/t ;                                        %reg. parameter decrease
    g = z - gamma*A*z ;                                       %gradient step
    uu = Projdualhinge(lambda*g,gamma*lambda)/lambda;         % proximal step ;   %
    z= uu + (t-alpha/4-1)*(uu-u)/(t+3*alpha/4) ;                            % inertial iterate
    b= toc ;
    %D(t) = (norm(X'*uu,2)^2)/(2*lambda^2) + sum(uu)/lambda;    %update dual objective
    %d(t) = norm(uu-u,2)^2 ;                                    % finite difference
    U(:,t)=uu ;                                        %update dual 
    ww= -X'*uu ; W(:,t)= ww;                           % update primal iterate, matrix
    M(t)=min(X*ww)/norm(ww);  wn(:,t)=ww/norm(ww);     % update normalized margin, error from wnast
    time(t-1)=b ;
    u = uu; 
end

end