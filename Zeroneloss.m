function [L] = Zeroneloss(K,u)
n=length(K(:,1));
L=sum((K*u<1))/n ;

end