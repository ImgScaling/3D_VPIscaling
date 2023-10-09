%  Test3 of the paper:  upscaling for 2 the Dicom image in the directory
% 'EsempioDICOM'
% In this test, two directions remain unchanged, and we upscale only along
% the remaining direction.
% To do this, we produce the Input  image, by  dropping slides with step d, on the chosen 
% direction and then we upscale with all the methods proposed in the paper
% Parameter to be set:%
% direction: fix the direction along which we upscale
% slice: the number of the slice you wanto to graphic
clear
A=dicomCollection('EsempioDICOM');D=dicomreadVolume(A);D = squeeze(D);
I_target=squeeze(D);
s=size(I_target);% s is the array of the size of the  target image axbxc to [a,b,c]
% in this test the sizes are [512 512 361]
%
% choose the direction x=1,y=2,z=3 and the slice you wanto to show
direction=3; slice=99;
% in other test you can choose the step d to use to decimate the image along the direction
% In this test d=2;
d=2;
switch direction
    case 1
    I_in=I_target(1:d:s(1),:,:);
    case 2
    I_in=I_target(:,1:d:s(2),:);
    case 3
        I_in=I_target(:,:,1:d:s(3));
    otherwise
        error(' check the direction to thin')   
end
%VPI supervised method.
%  m is seleceted as m=fix(n*theta), with theta \in (0,1)
k=1;
for theta=0.1:0.1:0.9
    [I_fin] = VPI_dicom(I_in,theta,s); I_fin=uint16(I_fin);
    ps_VPI(k)=psnr(I_target,I_fin);
    sim_VPI(k)=ssim(I_target,I_fin);
    vet_theta(k)=theta;
    k=k+1;
end
[ps_VPI_opt, ind]=max(ps_VPI);sim_VPI_opt=max(sim_VPI);
%
% BIC and Lanczos methods
I_finBIC=imresize3(I_in,s);I_finLanc=imresize3(I_in,s,'lanczos3');
ps_BIC=psnr(I_target,I_finBIC);ps_lanc=psnr(I_target,I_finLanc);
sim_lanc=ssim(I_target,I_finLanc);sim_BIC=ssim(I_target,I_finBIC);
%
% LCI method
[I_fin_lag] = VPI_dicom(I_in,0,s); I_fin_lag=uint16(I_fin_lag);
ps_lag=psnr(I_target,I_fin_lag); sim_lag=ssim(I_target,I_fin_lag);
% visualization of 3D results
disp('  3D results:      PSNR and   SSIM  ')
disp('BIC     Lanczos      Lagrange      VPI')
disp([ps_BIC,ps_lanc,ps_lag, ps_VPI_opt]);
disp([sim_BIC,sim_lanc,sim_lag,sim_VPI_opt]);
theta=vet_theta(ind);% the variable theta contains the optimal theta in (0,1) s.t. 
% m=fix(n*theta)
%
% This second part is devoted to obtain 2d PSNR and SSIM, for each 2d slide
 [I_fin] = VPI_dicom(I_in,theta,s); I_fin=uint16(I_fin);
 j=1;
 for k=2:d:s(direction) 
    switch direction
       case 1
       C1= I_target(k,:,:); C2=I_fin(k,:,:);C3=I_finBIC(k,:,:);
       C4=I_fin_lag(k,:,:);C5=I_finLanc(k,:,:);
         case 2
       C1= I_target(:,k,:); C2=I_fin(:,k,:);C3=I_finBIC(:,k,:);
       C4=I_fin_lag(:,k,:);C5=I_finLanc(:,k,:);
         case 3
       C1= I_target(:,:,k); C2=I_fin(:,:,k);C3=I_finBIC(:,:,k);
       C4=I_fin_lag(:,:,k);C5=I_finLanc(:,:,k);    
   end
 p_vp(j)=psnr(C1,C2);sim_vp(j)=ssim(C1,C2);
 p_bic(j)=psnr(C1,C3);sim_bic(j)=ssim(C1,C3);
 p_lag(j)=psnr(C1,C4); sim_lag(j)=ssim(C1,C4);
 p_lanc(j)=psnr(C1,C5);  sim_lanc(j)=ssim(C1,C5);
 j=j+1;
 end
disp('mean PSNR      and           mean SSIM in 2D ')
disp('       BIC   Lanczos     Lag    VPI  ')
matricep=[p_bic' p_lanc' p_lag' p_vp'];
matrices=[sim_bic' sim_lanc' sim_lag' sim_vp'];
Q=[mean(p_bic), mean(p_lanc), mean(p_lag),mean(p_vp)];
T=[mean(sim_bic),  mean(sim_lanc), mean(sim_lag),mean(sim_vp)];
disp(Q)
disp(T)
% graphic of the specific slices
% set the number of the slice to show
dir=[0 0 0];dir(direction)=1;
sliceViewer(I_target,'SliceDirection',dir, 'SliceNumber',slice)
figure
sliceViewer(I_fin,'SliceDirection',dir, 'SliceNumber',slice) % VPI
figure
sliceViewer(I_finBIC,'SliceDirection',dir, 'SliceNumber',slice)% BIC
sliceViewer(I_finLanc,'SliceDirection',dir, 'SliceNumber',slice)% Lanczos
sliceViewer(I_fin_lag,'SliceDirection',dir, 'SliceNumber',slice) % LCI
% 3D visualization 
volumeViewer(D)




  