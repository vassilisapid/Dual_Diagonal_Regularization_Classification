close all
clear
clc

%oldpath = addpath(fullfile(matlabroot,'examples','nnet','main'));
filenameImagesTrain = 'train-images-idx3-ubyte';
filenameLabelsTrain = 'train-labels-idx1-ubyte';
filenameImagesTest = 't10k-images-idx3-ubyte';
filenameLabelsTest = 't10k-labels-idx1-ubyte';

XTrain = processMNISTimages(filenameImagesTrain);
YTrain = processMNISTlabels(filenameLabelsTrain);
XTest = processMNISTimages(filenameImagesTest);
YTest = processMNISTlabels(filenameLabelsTest);
Xtrain=reshape(XTrain,[28,28,60000]);
Ztrain=reshape(XTrain,[784,60000])';
dd=28; 

%XTrain01 = processMNISTimages(filenameImagesTrain);
idtr0=find(YTrain=="3"); Ytr0=YTrain(idtr0); idtr1=find(YTrain=="5"); Ytr1=YTrain(idtr1); Ytr01=[Ytr0;Ytr1];
Xtr0=XTrain(:,:,idtr0); Xtr1=XTrain(:,:,idtr1); 
Xtr01=zeros(dd,dd,length(Ytr01)); Xtr01(:,:,[1:length(Ytr0)])=Xtr0; Xtr01(:,:,[length(Ytr0)+1:end])=Xtr1;
ntr01=length(Ytr01);

idts0=find(YTest=="3"); Yts0=YTest(idts0); idts1=find(YTest=="5"); Yts1=YTest(idts1); Yts01=[Yts0;Yts1]; nts01=length(Yts01);
Xts0=XTest(:,:,idts0); Xts1=XTest(:,:,idts1);
Xts01=zeros(dd,dd,length(Yts01)); Xts01(:,:,[1:length(Yts0)])=Xts0; Xts01(:,:,[length(Yts0)+1:end])=Xts1;
Ztr=reshape(Xtr01,[784,ntr01])'; Zts=reshape(Xts01,[784,nts01])';

[o,~,Gtr]=unique(Ytr01); [o,~,Gts]=unique(Yts01); Gtr(Gtr==1)=0; Gtr(Gtr==2)=1; Gts(Gts==1)=0; Gts(Gts==2)=1;