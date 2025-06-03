function [uu,M,Lts,time]=GDexpvarstepKernel(u0,Ats,A,gamma,T,Yts)
 
%Kts= KernelMatrix(Xts, X, kernel, kerpar); B=-Yts.*(Kts.*Y); 
Lts=ones(T,1)*calcErr(sign(-Ats*u0),Yts);
M=min(-A*u0)*ones(T,1);  
time=zeros(T-1,1) ;u=u0 ;
for t=2:T
    tic ;
    v=exp(-A*u);
    g = -exp(-A*u)  ;                  %gradient of Logistic                  
    uu = u - gamma*g /(sqrt((t)*dot(A*v,v))) ;
    b=toc ;
    ypred=sign(Ats*uu); 
    Lts(t)=calcErr(ypred,Yts);
    M(t)=min(A*uu)/sqrt(dot(A*uu,uu));   
    time(t-1)=b ;
    u = uu; 
end

end