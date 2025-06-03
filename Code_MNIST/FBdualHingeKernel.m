function [u,M,Lts,time]=FBdualHingeKernel(u0,Ats,A,gamma,lambda0,T,Yts)
%[d,n]=size(X');
% Initialization8

Lts=ones(T,1)*calcErr(sign(-Ats*u0),Yts);
M=min(-A*u0)*ones(T,1);                                               %Margin                 
time=zeros(T-1,1) ;
u=u0;
for t=2:T
    tic ;
    lambda=lambda0/(t) ;                                                     % reg. parameter decrease
    %gg= u -gamma*(A*u +ones(n,1));
    %uu=Projectionlambda(gg,lambda);
    g = u - gamma*A*u ;                                           % gradient step  
    uu = Projdualhinge(lambda*g,gamma*lambda)/lambda;                      % Projection step
    b=toc ; 
    ypred=sign(-Ats*uu);
    Lts(t)=calcErr(ypred,Yts);
    %Lts(t)=Zeroneloss(B,uu) ; 
    %Ltr(t)=Zeroneloss(-A,uu) ;
    M(t)=min(-A*uu)/sqrt(dot(A*uu,uu));                                                        %Margin update      
    time(t-1)=b ;
    u = uu; 
end

end