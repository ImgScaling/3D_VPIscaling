% script Test 3
clear
s=1/2;fatt=1/s;jj=1;
n=128;
I_target=phantom3d('Shepp-Logan',n);I_target=uint8(I_target);
I_in=phantom3d('Shepp-Logan',fix(n*fatt));I_in=uint8(I_in);
I_finMATLAB=imresize3(I_in,[n n n],'linear');I_fin_bic=uint8(I_finMATLAB);
I_finLanc=imresize3(I_in,[n n n],'lanczos3');I_fin_lanc=uint8(I_finLanc);
I_fin_lag = VPI_dicom(double(I_in),0,[n n n]);I_fin_lag =uint8(I_fin_lag);
I_target=uint8(I_target);
  k=1;
for theta=.05:.05:0.9
    I_fin = VPI_dicom(double(I_in),theta,[n n n]);I_fin=uint8(I_fin);
    ps_nostro(k)=psnr(I_target,I_fin);
    k=k+1;
end
ps_vpi=max(ps_nostro);
ps_bic=psnr(I_target,I_fin_bic);sqrt(immse(I_target,I_fin_bic))
ps_lanc=psnr(I_target,I_fin_lanc);
ps_lag=psnr(I_target,I_fin_lag);sqrt(immse(I_target,I_fin_lag))
ris(1,1:4)=[ps_bic, ps_lanc, ps_lag,ps_vpi];
ris1=['     Bic  ',    '   Lanczos     Lag         VP  '];
disp(ris1)
disp(ris)
