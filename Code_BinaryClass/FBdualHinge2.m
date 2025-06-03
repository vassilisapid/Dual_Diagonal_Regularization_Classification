function [u,U,ww,W,wn,M,time]=FBdualHinge2(w0,X,gamma,lambda0,T)
d=length(w0); n=length(X(:,1)); 

% Initialization
u=X*w0 ; U=u.*ones(n,T); W=w0.*ones(d,T);               % u0 dual, dual matrix, primal matrix
wn=w0.*ones(d,T);                                       % normalized error from solution wnast

A=X*(X');                          %Gradient

%D = ( (norm(-X'*u)^2)/2 + sum(lambda0*u)/lambda0 ).*ones(T,1) ;      %Dual objective
M=min(X*w0).*ones(T,1);                                               %Margin                 
time=zeros(T-1,1) ;

for t=2:T
    tic ;
    lambda=lambda0/t ;                                             % reg. parameter decrease
    g = u - gamma*A*u ;                                            % gradient step  
    uu = Projdualhinge(lambda*g,gamma*lambda)/lambda;              % Projection step
    b=toc ;
    %d(t) = norm(uu-u,2)^2 ;                                       % finite difference
    %D(t) = (norm(X'*uu,2)^2)/(2*lambda^2) + sum(uu)/lambda;       % update Dual objective
    
    U(:,t)=uu ;                                 %Dual update matrix  
    ww= -X'*uu ;                                %Primal update
    W(:,t)= ww;                                 %Primal update matrix,
    wn(:,t)=ww/norm(ww);                        %prime normalized update
    M(t)=min(X*ww)/norm(ww);                    %Margin update
    %b(t)=-(min(X(1:n/2,:)*ww)+max(X(n/2:end,:)*ww))/2 ;       
    
    time(t-1)=b ;
    u = uu; 
end

end