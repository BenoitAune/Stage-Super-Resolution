% % Super-resolution demo script as explained in: 
% 
% S. Parameswaran, C-A. Deledalle, L. Denis and T. Q. Nguyen, "Accelerating
% GMM-based patch priors for image restoration: Three ingredients for a 
% 100x speed-up", arXiv.
%
% We use EPLL and FEPLL
% I compile with differents K numbers of components and differents q
% magnification factors.
% I compile 10 times to have stable results and for 5 images with size
% 512*512


clear all
close all
addpathrec('.')
deterministic('on');

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
% Parameters
sig    = 2;
op     = operators('subresolution', M, N, 'width', 0.5, 'factor', 1/q);
sig    = sig/255;
y      = op.A(x) + sig * randn(op.osize);

% Load prior computed offline
prior_model = get_prior_model();
nu=ones(prior_model.GS.dim,1)*2;
% Run FEPLL
tstart = tic;
xhat = fepll(y,q, sig, prior_model, 'operator', op);
tFEPLL = toc(tstart);
for k=1:prior_model.GS.nmodels
    bla{k}=prior_model.GS.mu(:,k);
    prior_model.GS.nu{k}=nu;
end
prior_model.GS.mu=bla;
prior_model.name='gmm';

tstart = tic;
xhat2 = ggmm_epll2(y,q, sig, prior_model, 'operator', op);
tEPLL = toc(tstart);

% Display
figure()
subplot(2,3,1)
imagesc(x, [0 1]);title('Image originale X');
subplot(2,3,2)
imagesc(y, [0 1]);title('Image observée Y (connue)');
subplot(2,3,4)
xb = imresize(y, op.isize, 'bicubic');
imagesc(xb, [0 1]);title(sprintf('Z: PSNR %.1f SSIM %.3f', ...
              psnr(xb, x), ...
              ssim(xb, x)));
subplot(2,3,5)
imagesc(xhat2, [0 1]);title(sprintf('X_h: (EPLL): PSNR %.1f SSIM %.3f time %.2f', ...
              psnr(xhat2, x), ...
              ssim(xhat2, x), tEPLL));
subplot(2,3,6)
imagesc(xhat, [0 1]);title(sprintf('X_h: (FEPLL): PSNR %.1f SSIM %.3f time %.2f', ...
              psnr(xhat, x), ...
              ssim(xhat, x), tFEPLL));
colormap gray
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 50 30]);
saveas(gcf,['fig/',img,'_fepll_q',num2str(q)],'fig');
saveas(gcf,['png/',img,'_fepll_q',num2str(q)],'png');
tab=[psnr(xhat2, x),ssim(xhat2, x) , tEPLL;psnr(xhat, x),ssim(xhat, x),tFEPLL];
tabq = cat(2,tabq,tab);
end
tabimg=cat(1,tabimg,tabq);
end
csvwrite(['tab/tab_fepll','.csv'],tabimg) ;

