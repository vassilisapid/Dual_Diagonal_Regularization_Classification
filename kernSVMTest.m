function ypred = kernSVMTest(c, Xtr,Ytr, kernel, kerpar, Xts)
    K = KernelMatrix(Xts, Xtr, kernel, kerpar);
    ypred = sign(-(K*diag(Ytr))*c);
end