%proximal of dual of a*hinge Prox_{al*}(u)  

function [p]=Projdualhinge(u,a)
%d=ones(length(w),1);

p= -1.*(u<=a-1) + (u-a).*(a-1<u & u<a) +0.*(u>=a) ;

end