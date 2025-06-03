function [p]=Projectionlambda(u,a)
%d=ones(length(w),1);

p= -(u<=-1/a)/a + (u).*(-1/a<u & u<0) +0.*(u>=0) ;

end