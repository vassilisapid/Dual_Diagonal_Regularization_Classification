function [u,M,Lts]=GDexpPrimalDual(u0,Ats,A,theta,T,Yts)
n=length(A(:,1));     
Lts=ones(T,1)*calcErr(sign(-Ats*u0),Yts);
M=min(-A*u0).*ones(T,1); G=0*ones(1,1);
beta=0 ; q=ones(n,1)/n; u = u0 ; 
for t=2:T
    %S=sum(exp(-u))/n; G= -A*exp(-u)/n ;                %gradient of Exp
    G= beta*(G+q) ;                               %inertial term
    uu = u - theta*(G + q) ;                      %update
    q=exp(A*uu)/(norm(exp(A*uu),1));
    ypred=sign(-Ats*uu);
    Lts(t)=calcErr(ypred,Yts);
    M(t)=min(-A*uu)/sqrt(dot(A*uu,uu));
    u = uu;  
end

end