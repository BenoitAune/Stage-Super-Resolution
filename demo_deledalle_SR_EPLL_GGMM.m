% % Denoising demo script as explained in:
%
% C-A. Deledalle, S. Parameswaran, and T. Q. Nguyen, "Image
% restoration with generalized Gaussian mixture model patch
% priors", arXiv.
%
% We use EPLL with GMM, LMM and GGMM
% I compile with differents K numbers of components and differents q
% magnification factors.
% I compile 10 times to have stable results and for 5 images with size
% 512*512
cd(fileparts(mfilename('fullpath')));
addpathrec('.')
deterministic('on');

clear all
close all
cd(fileparts(mfilename('fullpath')));
addpathrec('.')
deterministic('on');
% Parameters

Q=[2,3];


tabimg=[];
for im=1:5
    if im==1
    img='lena';
    ext='.tif';
    x = rgb2gray(double(imread(['images/',img,ext]))/255);
    else
        if im==2
             img='camera';
            ext='.jpg';
            x = rgb2gray(double(imread(['images/',img,ext]))/255);
        else
            if im==3
            img='barbara';
            ext='.png'; 
            x =double(imread(['images/',img,ext]))/255;
            else
                if im==4
                img='hill';
                ext='.png';
                x =double(imread(['images/',img,ext]))/255;
                else
                img='shape';
                x=ones(512,512)*1;
                x(50:150,50:150)=0;
                x(350:450,350:450)=0;
                for s=0:0.01:2*pi
                   x(-round(150/2*cos(s))+100:round(150/2*cos(s))+100,-round(150/2*sin(s))+400:round(150/2*sin(s))+400)=0;
                end
                for s=50:50:150
                   x(300:500,s-12:s+12)=0; 
                   x(s-12+300:s+12+300,20:200)=0;
                end   
                end
            end
        end
    end    
[M, N] = size(x);
tabq=[];
for q=Q
sig    = 2;  
% Load and generate images
op     = operators('subresolution', M, N, 'width', 0.5, 'factor', 1/q);
sig    = sig / 255;
y      = op.A(x) + sig * randn(op.osize);

% Load prior computed offline
prior_model{1} = get_prior('gmm2');
prior_model{2} = get_prior('lmm');
prior_model{3} = get_prior('ggmm');

% Run GGMM EPLL
for k = 1:length(prior_model)
    tstart = tic;
    prior_model2=prior_model{k};
    xhat{k} = ggmm_epll2(y,q, sig, prior_model2, 'operator', op);
    toc(tstart);
end
prior_model{1}.name='gmm';

% Display
figure()
colormap gray
subplot(2,3,1)
imagesc(x, [0 1]);title('Original image X');
subplot(2,3,2)
imagesc(y, [0 1]);title(' Observed image Y (known)');
subplot(2,3,3)
xb = imresize(y, op.isize, 'bicubic');
imagesc(xb, [0 1]);title(sprintf('Z: (PSNR %.1f SSIM %.3f)', ...
              psnr(xb, x), ...
              ssim(xb, x)));
tab=[];
for k = 1:length(prior_model)
    subplot(2,3,3+k)
    imagesc(xhat{k}, [0 1]);
    title(sprintf('X_h: %s+EPLL (PSNR %.2f, SSIM %.3f)', ...
                  upper(prior_model{k}.name), ...
                  psnr(xhat{k}, x), ...
                  ssim(xhat{k}, x)));
    tabk=[psnr(xhat{k},x),ssim(xhat{k}, x)];  
    tab=cat(1,tab,tabk);
end
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 50 30]);
saveas(gcf,['fig/',img,'_ggmm_epll_q',num2str(q)],'fig');
saveas(gcf,['png/',img,'_ggmm_epll_q',num2str(q)],'png');
tabq=cat(2,tabq,tab);
end
tabimg=cat(1,tabimg,tabq);
end
csvwrite(['tab/tab_ggmm_epll','.csv'],tabimg) ;
