close all
clear
clc

T=readtable('HTRU_2.csv');
Y=T(:,end); Y=table2array(Y); Y(Y==0)=-1;
X=T(:,1:end-1);  X=table2array(X); [n,d]=size(X);
M=mean(X); Sigma=std(X); Z=(X-M.*ones(n,d))./Sigma;

ptrain=0.7; Index = randperm(n,ceil(n*ptrain));
Xtr=X(Index,:); Ytr=Y(Index); Xts=X(setdiff(1:end,Index),:); Yts=Y(setdiff(1:end,Index));
