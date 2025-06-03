function [u,M,Lts]=FBexp(u0,Ats,A,gamma,lambda0,T,Yts)
%[d,n]=size(X');
% Initialization8

Lts=ones(T,1)*calcErr(sign(-Ats*u0),Yts);
M=min(-A*u0)*ones(T,1);                                               %Margin                 
u=u0;
for t=2:T
    lambda=lambda0/(2^(min(t,1021))) ;                                                     % reg. parameter decrease
    %gg= u -gamma*(A*u +ones(n,1));
    %uu=Projectionlambda(gg,lambda);
    g = u - gamma*A*u ;                                           % gradient step  
    uu = Projdualhinge(lambda*g,gamma*lambda)/lambda;                      % Projection step
    ypred=sign(-Ats*uu);
    Lts(t)=calcErr(ypred,Yts);
    %Ltr(t)=Zeroneloss(-A,uu) ;
    M(t)=min(-A*uu)/sqrt(dot(A*uu,uu));                                                        %Margin update      
    u = uu; 
end