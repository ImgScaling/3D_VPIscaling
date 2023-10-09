% imagine BRAIN, test 2 of the papers
clear
load MRI
D = squeeze(D);%I_in=D;
I_target=squeeze(D);
% s contains the size of the  target image [128, 128, 27]
s=[128,128,27];
% 
% ss has to be set and contains the size of the initial image
% the possible sizes in upscaling are 
% 
ss=[64 64 14];%ss=[43,43,9]; % ss=[32 32 7];
% the posible sizes in downscaling are
% ss=[256 256 54];%ss=[384 384 81];%ss=[512   512   108];
%
% Set the method  to be used to produce the initial image. 
choice=1;
switch  choice
    case  1
    I_in=imresize3(I_target,ss);% Bicubic
    case 2
    I_in=imresize3(I_target,ss,'lanczos3');% Lanczos
    case 3
    I_in=VPI_dicom(I_target,0,ss);I_in=uint8(I_in);% Lagrange
    case 4
    I_in=VPI_dicom(I_target,0.7,ss);I_in=uint8(I_in);% VPI
    otherwise
        error('case unknown')
end
% resize the initial image with the comparision methods BIC, Lanczos, LCI
% and evaluate psnr and ssim
I_fin_BIC=imresize3(I_in,s);
psnr_BIC=psnr(I_fin_BIC,I_target);sim_BIC=ssim(I_fin_BIC,I_target);

I_fin_lanc=imresize3(I_in,s,'lanczos3');
psnr_lanc=psnr(I_fin_lanc,I_target);sim_lanc=ssim(I_fin_lanc,I_target);


I_fin_lag=VPI_dicom(I_in,0,s);I_fin_lag=uint8(I_fin_lag);
psnr_lag=psnr(I_fin_lag,I_target);sim_lag=ssim(I_fin_lag,I_target);
% show the obtained results by the comparison methods
disp(' BIC          Lanczos       LCI')
disp([psnr_BIC,  psnr_lanc,   psnr_lag])
disp([sim_BIC,  sim_lanc,   sim_lag])
% supervised VPI method selecting theta in (0,1), s.t. m=fix(n*theta)
% 
k=1;theta_vet=0.15:0.1:0.9;
for theta=theta_vet
[I_fin] = VPI_dicom(I_in,theta,s);I_fin=uint8(I_fin);
psnr_VPI(k)=psnr(I_fin,I_target);
sim_VPI(k)=ssim(I_fin,I_target);
k=k+1;
end
[psnr_VPI_opt,ind]=max(psnr_VPI); sim_VPI_opt=max(sim_VPI);
disp( '    VPI results   ')
disp ( [psnr_VPI_opt, sim_VPI_opt])
theta_opt=theta_vet(ind);mm=fix(theta_opt*ss);
disp('         optimal m')
disp(mm)
