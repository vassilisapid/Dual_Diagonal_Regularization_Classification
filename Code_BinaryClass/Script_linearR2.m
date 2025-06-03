%Classification in R^2

clc;
clear all;
close all;
N=40;                   %Number of data-points in each class
xsupp1=[1/2,3/2]; xsupp2 = [3/2,1/2]; xsupp3=-[1/2,3/2]; xsupp4 = -[3/2,1/2];       %xsupp : support vectors 

[X, Y] = MixGauss([[2;2], [-2;-2]],[0.45, 0.45], N);                          % random data generation
X(1,:)=xsupp1 ; X(2,:)=xsupp2 ; X(2*N-1,:)=xsupp3 ; X(2*N,:)=xsupp4 ;
Y=2*mod(Y, 2)-1;  
wast=[1;1]/2 ; Mast=min(Y.*(X*wast)); wnast=wast/norm(wast,2);      % wast : min norm solution ,  Mast : max-mmargin



A=Y.*X ; L=max(eig(A*(A'))) ; mu=min(abs(eig(A*(A'))));

w0=[0;0]; gamma=1/L;  T=1000; lambda0=0.001;    % w0 : Start. point, gamma : steps size, lambda : reg parameter

[u1,U1,w1,W1,wn1,M1,t1]=FBdualHinge2(w0,A,1.4*gamma,lambda0,T);    %FB on the dual

alpha1=10 ;
[u2,U2,w2,W2,wn2,M2,t2]=iFBdualHinge2(w0,A,alpha1,gamma,lambda0,T);   %i-FB on the dual with alpha 

alpha2=30 ;
[u3,U3,w3,W3,wn3,M3,t3]=iFBdualHinge2(w0,A,alpha2,gamma,lambda0,T);

alpha3=50 ;
[u4,U4,w4,W4,wn4,M4,t4]=iFBdualHinge2(w0,A,alpha3,gamma,lambda0,T);

q1=zeros(T,1); q2=q1 ; q3=q1 ; q4=q1 ;
for t=1:T
    q1(t)=1-dot(W1(:,t)/norm(W1(:,t)),wnast);
    q2(t)=1-dot(W2(:,t)/norm(W2(:,t)),wnast);
    q3(t)=1-dot(W3(:,t)/norm(W3(:,t)),wnast);
    q4(t)=1-dot(W4(:,t)/norm(W4(:,t)),wnast);
end    

t = linspace(-3,3,100)'; 

figure(1)
scatter(X(1:N,1),X(1:N,2),40,Y(1:N),'filled') ; hold on
scatter(X(N:2*N,1),X(N:2*N,2),40,Y(N:2*N),'filled') 
%scatter3(X(:,1),X(:,2),X(:,3),40,Y,'filled') ; title('dataset1'); hold on;
%plot(t,-t+2,'.k','LineWidth',0.8);  hold on;
%plot(t,-t-2,'.k','LineWidth',0.8);
plot(t,t,'k','LineWidth',1); hold on
%fimplicit(f1,'--g','LineWidth',2) ; 
%fimplicit(f2,'--r','LineWidth',1.5) ;  fimplicit(f3,'--m','LineWidth',1.5) ; 
%fimplicit(fvNexp,'--y','LineWidth',1) ; fimplicit(fexp,'--b','LineWidth',1) ;  
plot(t,-(w1(1,end))/w1(2,end)*t,'g-o','MarkerIndices',1:9:length(t),'LineWidth',1.2) ; 
plot(t,-w2(1,end)/w2(2,end)*t,'r-+','MarkerIndices',1:7:length(t),'LineWidth',1.2) ;
plot(t,-w3(1,end)/w3(2,end)*t,'m-*','MarkerIndices',1:5:length(t),'LineWidth',1.2) ;  
plot(t,-w4(1,end)/w4(2,end)*t,'b-.','MarkerIndices',1:13:length(t),'LineWidth',1.2) ;  
%plot(t,-wvNexp(1,end)/wvNexp(2,end),'y--','LineWidth',1) ;
%plot(t,-wexp(1,end)/wexp(2,end)*t,'b--','LineWidth',1) ;
xlabel('x')
ylabel('$$ y$$','Interpreter','Latex')
legend('data$$+$$','data$$-$$','$$w_*$$','FB','iFB $(\alpha=10)$','iFB $(\alpha=30)$','iFB $(\alpha=50)$','Interpreter','Latex')
axis equal
grid on



figure(3)
plot(log(vecnorm(wn1-wnast,1)),'g','LineWidth',2) ; hold on;
plot(log(vecnorm(wn2 -wnast,1)),'r','LineWidth',2)
plot(log(vecnorm(wn3-wnast,1)),'m','LineWidth',2)
plot(log(vecnorm(wn4-wnast,1)),'b','LineWidth',2)
%loglog(wnvNexp,'y','LineWidth',1.2)
%loglog(wnsub,'b','LineWidth',1)
%loglog(wnexp,'k','LineWidth',1)
%loglog(wn6,'c','LineWidth',1)
title('Error values') 
xlabel('Iterations')
ylabel('$$ \ln\|\bar{w_t}-\bar{w^\ast}\|$$','Interpreter','Latex')
legend('FB','iFB $(\alpha=10)$','iFB $(\alpha=30)$','iFB $(\alpha=50)$','Interpreter','Latex')
grid on

figure
plot(log(abs(M1-1/norm(wast))),'g','LineWidth',2) ; hold on;
plot(log(abs(M2-1/norm(wast))),'r','LineWidth',2);
plot(log(abs(M3-1/norm(wast))),'m','LineWidth',2);
plot(log(abs(M4-1/norm(wast))),'b','LineWidth',2);
title('Margin gap') 
xlabel('Iterations')
ylabel('$$ \mathbf{\ln}(\textbf{Margin gap})$$','Interpreter','Latex')
legend('FB','iFB $(\alpha=10)$','iFB $(\alpha=30)$','iFB $(\alpha=50)$','Interpreter','Latex')
grid on
figure('Name','Angle gap');
plot(log(q1),'g','LineWidth',2) ; hold on;
plot(log(q2),'r','LineWidth',2)
plot(log(q3),'m','LineWidth',2)
plot(log(q4),'b','LineWidth',2)
%plot(log(abs(MvNexp-Mast)),'y','LineWidth',1.2)
%plot(log(abs(Msub-Mast)),'b','LineWidth',1)
%plot(log(abs(Mexp-Mast)),'k','LineWidth',1)
%plot(log(abs(M6-Mast/norm(wast))),'c','LineWidth',1)
title('Angle gap') 
xlabel('Iterations')
ylabel('$$ \mathbf{\ln}(\textbf{Angle gap})$$','Interpreter','Latex')
legend('FB','iFB $(\alpha=10)$','iFB $(\alpha=30)$','iFB $(\alpha=50)$','Interpreter','Latex')
grid on