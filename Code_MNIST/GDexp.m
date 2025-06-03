function [w,W,d,wn,M,b]=GDexp(w0,X,wnast,gamma,T)
n=length(X(:,1)); w = w0 ; b=zeros(T,1);
W = (sum(exp(-X*w0))/log(2)).*ones(T,1) ; M=min(X*(w0)).*ones(T,1);
d=ones(T,1);  wn=norm(w0/norm(w0) - wnast,2).*ones(T,1);
for t=2:T
    g = -X'*exp(-X*w)   ;                  %gradient of Logistic                  
    ww = w - gamma*g  ; %/(sqrt(t+1)*norm(g)) ;
    W(t) = sum(exp(-X*w)) ;     %-sum(log(1+exp(-X*([1;1]/2))))/log(2) ; 
    d(t) = norm(ww-w,2)^2 ;
    w = ww; 
    wn(t)=norm( w/norm(w) - wnast,2); M(t)=min(X*w)/norm(w); b(t)=-(min(X(1:n/2,:)*w)+max(X(n/2:end,:)*w))/2 ;
end

end