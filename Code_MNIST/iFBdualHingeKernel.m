function [u,M,Lts,time]=iFBdualHingeKernel(u0,Ats,A,gamma,a,lambda0,T,Yts)
%[d,n]=size(X');
% Initialization  
Lts=ones(T,1)*calcErr(sign(-Ats*u0),Yts);
M=min(-A*u0)*ones(T,1);                                                %Margin                 
time=zeros(T-1,1) ;
u=u0; z=u0;
for t=2:T
    tic ;
    lambda=lambda0/(t) ;                                                     % reg. parameter decrease
    %gg= z -gamma*(A*z +ones(n,1));
    %uu=Projectionlambda(gg,lambda);
    g = z - gamma*A*z ;                                           % gradient step  
    uu = Projdualhinge(lambda*g,gamma*lambda)/lambda;                      % Projection step
    z=uu+(t-1)*(uu-u)/(t+a);
    b=toc ;
    ypred=sign(-Ats*uu);
    Lts(t)=calcErr(ypred,Yts);
    M(t)=min(-A*uu)/sqrt(dot(A*uu,uu));                                                        %Margin update      
    time(t-1)=b ;
    u = uu; 
end

end