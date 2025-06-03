function [u,M,Lts]=iGDexpPrimalDual(u0,Ats,A,theta,T,Yts)
n=length(A(:,1)); 

%Kts= KernelMatrix(Xts, X, kernel, kerpar); B=-Yts.*(Kts.*Y);
Lts=ones(T,1)*calcErr(sign(-Ats*u0),Yts);
 M=min(-A*u0).*ones(T,1); g=zeros(1,1);
q=ones(n,1)/n; u = u0 ; 
for t=2:T
    beta=t/(t+1) ;
    %S=sum(exp(A*u))/n; G= -A*exp(A*u)/n ;                %gradient of Exp
    g= beta*(g+q) ;                               %inertial term
    uu = u - theta*(g + q) ;                      %update
    %uu = u - theta*(g + G/S) ;
    q=exp(A*uu)/(norm(exp(A*uu),1));
    ypred=sign(-Ats*uu);
    Lts(t)=calcErr(ypred,Yts);
    M(t)=min(-A*uu)/sqrt(dot(A*uu,uu));
    u = uu;  
end

end