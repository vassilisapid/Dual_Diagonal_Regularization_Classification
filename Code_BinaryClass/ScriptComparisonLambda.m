%Classification in R^2 with gaussian kernel

clc;
clear all;
close all;

%% Data generation

ntr=300;                   %Number of data-points in each class
nts =300;
[Xtr, Ytr] = MixGauss([[-0.5;-0.5],[0.5;0.5]],[0.4,0.4],ntr);
[Xts, Yts] = MixGauss([[-0.5;-0.5],[0.5;0.5]],[0.4,0.4],nts);
Ytr(Ytr==2) = -1;
Yts(Yts==2) = -1;
% Adding Noise
p=0;                            % Noise level
Ytr = flipLabels(Ytr, p);       % Flipped labels
Yts = flipLabels(Yts, 0);
[n,d]=size(Xtr);

% Plotting
figure(1)
subplot(1,2,1)
scatter(Xtr(:,1),Xtr(:,2),30,Ytr,'filled');
title('Training Set');
subplot(1,2,2)
scatter(Xts(:,1),Xts(:,2),30,Yts,'filled');
title('Test Set');

%% Parameter selection
kernel='gaussian';
kerpar=0.12 ;
A= diag(Ytr)*KernelMatrix(Xtr, Xtr, kernel, kerpar)*diag(Ytr);
Kts=KernelMatrix(Xts, Xtr, kernel, kerpar); Ats=Kts*diag(Ytr);
L=max(eig(A)) ; %mu=min(abs(eig(A)));

u0=zeros(n,1)/100; gamma=1/L;  T=2000;    % w0 : Start. point, gamma : steps size, lambda : reg parameter

%% Algorithms
lambda0=80; 
[uu1,M1,H1ts]=FBdualHinge(u0,Ats,A,gamma,lambda0,T,Yts);    %FB on the dual
lambda01=10; 
[uu11,M11,H11ts]=FBdualHinge(u0,Ats,A,gamma,lambda01,T,Yts); 
lambda02=1; 
[uu12,M12,H12ts]=FBdualHinge(u0,Ats,A,gamma,lambda02,T,Yts); 
lambda03=0.01; 
[uu13,M13,H13ts]=FBdualHinge(u0,Ats,A,gamma,lambda03,T,Yts); 

alpha1=3 ;
[uu2,M2,H2ts]=iFBdualHinge(u0,Ats,A,gamma,alpha1,lambda0,T,Yts);
[uu21,M21,H21ts]=iFBdualHinge(u0,Ats,A,gamma,alpha1,lambda01,T,Yts); 
[uu22,M22,H22ts]=iFBdualHinge(u0,Ats,A,gamma,alpha1,lambda02,T,Yts); 
[uu23,M23,H23ts]=iFBdualHinge(u0,Ats,A,gamma,alpha1,lambda03,T,Yts); %i-FB on the dual with alpha 


%q1=zeros(T,1); q2=q1 ; q3=q1 ; q4=q1 ;
%for t=1:T
%    q1(t)=1-dot(W1(:,t)/norm(W1(:,t)),wnast);
%    q2(t)=1-dot(W2(:,t)/norm(W2(:,t)),wnast);
%    q3(t)=1-dot(W3(:,t)/norm(W3(:,t)),wnast);
%    q4(t)=1-dot(W4(:,t)/norm(W4(:,t)),wnast);
%end    


figure(2)
        scatter(Xtr(:,1),Xtr(:,2),30,Ytr,'filled');
        hold on
        title('kernelot')
        separatingKernSVM(uu1, Xtr, kernel, kerpar, Xtr, Ytr)

figure(3)
        scatter(Xtr(:,1),Xtr(:,2),30,Ytr,'filled');
        hold on
        title('kernelot2')
        separatingKernSVM(uu2, Xtr, kernel, kerpar, Xtr, Ytr)
        %separatingKernSVM(uvs, Xtr, kernel, kerpar, Xtr, Ytr)
        %separatingKernSVM(uPD, Xtr, kernel, kerpar, Xtr, Ytr)
        %separatingKernSVM(uiPD, Xtr, kernel, kerpar, Xtr, Ytr)
%Mast=MiPD(T)*ones(T,1);

figure(4)
plot(M1,'-g*','MarkerIndices',1:95:T,'LineWidth',1.9) ; hold on;
plot(M11,'--r+','MarkerIndices',1:90:T,'LineWidth',1.9);
plot(M12,'-.b^','MarkerIndices',1:90:T,'LineWidth',1.9);
plot(M13,'-.mo','MarkerIndices',1:95:T,'LineWidth',1.9);
title('Margin gap') 
xlabel('Iterations')
ylabel('$$ \textbf{Margin}$$','Interpreter','Latex')
legend('$\lambda_{0}=100$','$\lambda_{0}=10$','$\lambda_{0}=1$','$\lambda_{0}=0.01$','Interpreter','Latex')
grid on

figure(5)
plot(M2,'-g*','MarkerIndices',1:95:T,'LineWidth',1.9) ; hold on;
plot(M21,'--r+','MarkerIndices',1:90:T,'LineWidth',1.9);
plot(M22,'-.b^','MarkerIndices',1:90:T,'LineWidth',1.9);
plot(M23,'-.mo','MarkerIndices',1:95:T,'LineWidth',1.9);
title('Margin gap') 
xlabel('Iterations')
ylabel('$$ \textbf{Margin}$$','Interpreter','Latex')
legend('$\lambda_{0}=100$','$\lambda_{0}=10$','$\lambda_{0}=1$','$\lambda_{0}=0.01$','Interpreter','Latex')
grid on

figure(6)
plot(H1ts,'-g*','MarkerIndices',1:95:T,'LineWidth',1.9) ; hold on;
plot(H11ts,'--r+','MarkerIndices',1:90:T,'LineWidth',1.9);
plot(H12ts,'-.b^','MarkerIndices',1:90:T,'LineWidth',1.9);
plot(H13ts,'-.mo','MarkerIndices',1:95:T,'LineWidth',1.9);
title('Margin gap') 
xlabel('Iterations')
ylabel('$$ \textbf{Test error}$$','Interpreter','Latex')
legend('$\lambda_{0}=100$','$\lambda_{0}=10$','$\lambda_{0}=1$','$\lambda_{0}=0.01$','Interpreter','Latex')
grid on

figure(7)
plot(H2ts,'-g*','MarkerIndices',1:90:T,'LineWidth',1.9) ; hold on;
plot(H21ts,'--r+','MarkerIndices',1:90:T,'LineWidth',1.9);
plot(H22ts,'-.b^','MarkerIndices',1:95:T,'LineWidth',1.9);
plot(H23ts,'-.mo','MarkerIndices',1:95:T,'LineWidth',1.9);
title('Margin gap') 
xlabel('Iterations')
ylabel('$$ \textbf{Test error}$$','Interpreter','Latex')
legend('$\lambda_{0}=100$','$\lambda_{0}=10$','$\lambda_{0}=1$','$\lambda_{0}=0.01$','Interpreter','Latex')
grid on