function [tab_entier,tab_parthg,tab_parthd,tab_partbg,tab_partbd]=main_sandeep_SR_EM_MMSE_EPLL(K,q,tau,x_orig,img,init)

%input: K number of components
%       q the mignification factor
%       tau*tau the size of the patches of Y
%       x_orig the original image
%       img the name of the image
%       init the name of the initialisation
%output: tab_entier the tabular of the results in the full image
%       tab_parthg the tabular of the results in the part hg (top-left) of the image
%       tab_parthd the tabular of the results in the part hd (top-right) of the image
%       tab_partbg the tabular of the results in the part bg (bottom-left) of the image
%       tab_partbd the tabular of the results in the part bd (bottom-right) of the image

close all
colormap gray;

%Y=SHX+p
%Y the LR image , X the HR image , p a white gaussian noise ,
%S a subsampling operator and H a blur operator 

%% Parameters 
sig=2;         % standard noise
sig=sig/255;
h=37;          % Cluster dimension
gamma = 0.0000001; % regularization parameter to patches recovery

tau2=floor(tau/2);
tau2q=floor(tau*q/2);
epsilon = 0.000001; % the regularization parameter for the covariance
%% Load and generate Images

[Mx_orig, Nx_orig] = size(x_orig); 
op = operators('subresolution', Mx_orig, Nx_orig, 'width', 0.5, 'factor', 1/q); % S & H operator
y_obs      = op.A(x_orig) + sig * randn([Mx_orig/q, Nx_orig/q]); % y_obs is the subsamping of the HR image x_orig (LR observed image)

options    = makeoptions('operator',op);
autonorm   = getoptions(options, 'autonorm',      true);
m = min(y_obs(:));
M = max(y_obs(:));
a = 5;
if m < 0 - a*sig || M > 1 + a*sig
    if autonorm
        warning(['Image was not in the range [0 1], ' ...
                 'it has been automatically normalized. ' ...
                 'Note that this may alter the performances. ' ...
                 'For optimal results, please normalized the input ' ...
                 'image correctly in the range [0 1].']);
        shift = m+a*sig;
        scale = (M-m - 2*a*sig);
        sig = sig / scale;
        y_obs     = (y_obs - shift) / scale;
    else
        warning(['Image is not in the range [0 1]. ' ...
                 'Note that this may alter the performances. ' ...
                 'For optimal results, please normalized the input ' ...
                 'image correctly in the range [0 1], or ' ...
                 'rerun GGMM_EPLL with option "autonorm: true"']);
    end
end

[My_obs, Ny_obs] = size(y_obs);
wM=My_obs/2;    
wN=Ny_obs/2; 

%  randM=84/q;
%  randN=150/q;
randM=1;
randN=1;
%y is the image x in LR
y = y_obs(randM:randM+wM-1, randN:wN+randN-1);
[My, Ny] = size(y);
wM2 = wM*q;
wN2 = wN*q; % the size of part of the image X that we know
%x is the 'random' selection of a part of the HR image x_orig (HR observed image)
randM2=randM*q-(q-1);  
randN2=randN*q-(q-1);
x = x_orig(randM2:randM2+wM2-1, randN2:wN2+randN2-1);
[Mx, Nx] = size(x);



%% Patches collections
masky_obs   = getmask(My_obs, Ny_obs, tau,1, false);
yi_obs = getpatches(y_obs, tau, masky_obs);
[Myi_obs,~]= size(yi_obs);
masky   = getmask(My, Ny, tau, 1, false);
yi = getpatches(y, tau, masky);
[Myi,~]=size(yi);
maskx   = getmask(Mx, Nx, tau*q, q, false);
xi = getpatches(x, tau*q, maskx);
[Mxi,Nxi]=size(xi); %Mxi the number of pixels in a patch & Nxi the number of patches

vi= [xi;yi]; %the concatenate vector with the same number of patches
Mvi=Mxi+Myi; 

%% GMM parameters initialisation
% GMM parameters initialisation 

wVk=ones(1,K)*(K/Nxi);
hpatches=randi(Nxi,1,K);
vhnn = knnsearch(vi',vi(:,hpatches)','K',h,'Distance','euclidean');
muVk=zeros(Mvi,K);          %muVK = [muHk';muLk'] for HR and LR parts respectively
sigmaVk=zeros(Mvi,Mvi,K);    %sigmaVk = [sigmaHk sigmaHLk ; sigmaHLk' sigmaLk]

% parameters initialisation
for k=1:K
muVk(:,k)=mean(vi(:,vhnn(k,:)),2);    
sigmaVk(:,:,k)=(cov(vi(:,vhnn(k,:))'))+eye(Mvi)*epsilon;
end
 

model.mu=muVk;
model.Sigma=sigmaVk;
model.w=wVk;

%we use EM on vi to estimate all the parameters
%% EM
tstartEM = tic;
[~,model,~]=mixGaussEm(vi,model); %EM
tEM = toc(tstartEM);
muVk=model.mu;
sigmaVk= model.Sigma;
wVk=model.w;

tstart = tic;
muHk= zeros(Mxi,K);
sigmaHk=zeros(Mxi,Mxi,K);
for k=1:K
    muHk(:,k)=muVk(1:Mxi,k);
    sigmaHk(:,:,k)=sigmaVk(1:Mxi,1:Mxi,k);
end

%% Restoration with MMSE with 1 component
1
[xi_hat,k_hat] = restoration_MMSE(yi_obs,wVk,muVk,sigmaVk,Mxi); % MMSE
maskx_orig=zeros(My_obs*q,Ny_obs*q);
if mod(tau*q,2)==0
    maskx_orig(tau2q+1:q:My_obs*q-tau2q+1,tau2q+1:q:Ny_obs*q-tau2q+1)=1;
else
    maskx_orig(tau2q+1:q:My_obs*q-tau2q,tau2q+1:q:Ny_obs*q-tau2q)=1;
end
% reconstruction
X_hat = get_image2D(xi_hat,tau*q,My_obs*q,Ny_obs*q,k_hat,muVk(1:Mvi-Myi_obs,:),sigmaVk(1:Mvi-Myi_obs,1:Mvi-Myi_obs,:),maskx_orig,gamma);

% maskx_orig  = getmask(Mx_orig, Nx_orig, tau*q, q, randommask);
% [Mxi_orig,Nxi_orig]=size(xi);
% X_hat = projpatches(xi_hat,Mx_orig,Nx_orig,maskx_orig); %get image

%time
tMMSE=tEM + toc(tstart);

%% Restoration with MMSE with all components
2
tstart=tic;
G=fspecial('gaussian',5,0.5);
[Sub,Subt]=sub_operator(q,tau);
H=blur_operator(G,tau*q);
var=diag(ones(1,tau*tau)*(sig^2));
[xi_hat,k_hat] = restoration_MMSE_GMM(yi_obs,wVk,muVk,sigmaVk,var,Sub,Subt,H,q,tau); % MMSE 
%reconstruction
X_hat2 = get_image2D(xi_hat,tau*q,My_obs*q,Ny_obs*q,k_hat,muVk(1:Mvi-Myi_obs,:),sigmaVk(1:Mvi-Myi_obs,1:Mvi-Myi_obs,:),maskx_orig,gamma);
tEPLL=tEM + toc(tstart); %time MMSEGMM 

%% Restoration with EPLL
% 2
% tstart=tic;
% prior_model.name='gmm';
% nu=ones(tau*tau*q*q,1)*2;
% for k=1:K
%    [U , S] = eig(sigmaHk(:,:,k));
%    prior_model.GS.U{k}=U;
%    prior_model.GS.S{k}=S(S~=0);
%    prior_model.GS.wts(k)=wVk(k);
%    prior_model.GS.mu{k}=muHk(:,k);
%    prior_model.GS.nu{k}=nu;
% end
% prior_model.GS.nmodels=K;
% prior_model.GS.dim=tau*tau*q*q;
% 
% X_hat2 = ggmm_epll(y_obs,maskx_orig,gamma,q,k_hat, sig, prior_model, 'operator', op); % EPLL & get image
% 
% tEPLL=tEM + toc(tstart);

%% Display & Tab
3
tstart=tic;
xb = imresize(y_obs, op.isize, 'bicubic');
txb=toc(tstart);
psnr_X_hat1=psnr(X_hat,x_orig)
ssim_X_hat1=ssim(X_hat,x_orig)
psnr_X_hat2=psnr(X_hat2,x_orig)
ssim_X_hat2=ssim(X_hat2,x_orig)
psnrxb=psnr(xb, x_orig)
ssimxb=ssim(xb, x_orig)

img=[img,'_K', num2str(K),'q',num2str(q),'t',num2str(tau)];
subplot(1,2,1);imagesc(x_orig,[0 1]);title('Image originale X');
subplot(1,2,2);imagesc(x,[0 1]);title('Partie de X (connue)'); 
colormap gray

figure(2)
subplot(2,3,1);imagesc(x_orig,[0,1]);title('Image originale X');
subplot(2,3,2);imagesc(y_obs,[0 1]);title('Image observée Y (connue)');
subplot(2,3,4);imagesc(xb,[0 1]);
title(sprintf('Z: PSNR %.1f SSIM %.3f', psnr(xb, x_orig), ssim(xb, x_orig)));
% subplot(2,3,5);imagesc(X_hat,[0 1]);title(sprintf('X_h (MMSE): PSNR %.1f SSIM %.3f',psnr(X_hat, x_orig), ssim(X_hat, x_orig)));
% subplot(2,3,6);imagesc(X_hat2,[0 1]);title(sprintf('X_h (EPLL): PSNR %.1f SSIM %.3f',psnr(X_hat2, x_orig), ssim(X_hat2, x_orig)));
subplot(2,3,5);imagesc(X_hat,[0 1]);title(sprintf('X_h (MMSE with one component): PSNR %.1f SSIM %.3f',psnr(X_hat, x_orig), ssim(X_hat, x_orig)));
subplot(2,3,6);imagesc(X_hat2,[0 1]);title(sprintf('X_h (MMSE with all components): PSNR %.1f SSIM %.3f',psnr(X_hat2, x_orig), ssim(X_hat2, x_orig)));
colormap gray
tab_entier=[psnr_X_hat1, ssim_X_hat1, tMMSE; psnr_X_hat2, ssim_X_hat2, tEPLL; psnrxb, ssimxb,txb];
set(2, 'PaperUnits', 'centimeters');
set(2, 'PaperPosition', [0 0 35 30]);
saveas(2,['fig/',img,'_large_',init],'fig');
saveas(2,['png/',img,'_large_',init],'png');

figure()
posM=1;
posN=1;
posM2=1;
posN2=1;
subplot(2,3,1);imagesc(x_orig(posM:posM+wM2-1, posN:wN2+posN-1),[0,1]);title('Partie de X');
subplot(2,3,2);imagesc(y_obs(posM2:posM2+wM-1, posN2:wN+posN2-1),[0 1]);
title('Partie de Y (connue)');
subplot(2,3,4);imagesc(xb(posM:posM+wM2-1, posN:wN2+posN-1),[0 1]);
title(sprintf('z: PSNR %.1f SSIM %.3f',psnr(xb(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(xb(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))));
% subplot(2,3,5);imagesc(X_hat(posM:posM+wM2-1, posN:wN2+posN-1),[0 1]);
% title(sprintf('x_h (MMSE): PSNR %.1f SSIM %.3f',psnr(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))));
% subplot(2,3,6);imagesc(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1),[0 1]);
% title(sprintf('x_h (EPLL): PSNR %.1f SSIM %.3f',psnr(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))));
subplot(2,3,5);imagesc(X_hat(posM:posM+wM2-1, posN:wN2+posN-1),[0 1]);
title(sprintf('x_h (MMSE with one component): PSNR %.1f SSIM %.3f',psnr(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))));
subplot(2,3,6);imagesc(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1),[0 1]);
title(sprintf('x_h (MMSE with all components): PSNR %.1f SSIM %.3f',psnr(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))));
colormap gray
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 35 30]);
saveas(gcf,['fig/',img,'_parthg_',init],'fig'); % top-left
saveas(gcf,['png/',img,'_parthg_',init],'png');
tab_parthg=[psnr(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1));...
psnr(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1));...
psnr(xb(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(xb(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))];


figure()
posM=1;
posN=256;
posM2=1;
posN2=posN/q;
subplot(2,3,1);imagesc(x_orig(posM:posM+wM2-1, posN:wN2+posN-1),[0,1]);title('Partie de X');
subplot(2,3,2);imagesc(y_obs(posM2:posM2+wM-1, posN2:wN+posN2-1),[0 1]);
title('Partie de Y (connue)');
subplot(2,3,4);imagesc(xb(posM:posM+wM2-1, posN:wN2+posN-1),[0 1]);
title(sprintf('z: PSNR %.1f SSIM %.3f',psnr(xb(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(xb(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))));
% subplot(2,3,5);imagesc(X_hat(posM:posM+wM2-1, posN:wN2+posN-1),[0 1]);
% title(sprintf('x_h (MMSE): PSNR %.1f SSIM %.3f',psnr(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))));
% subplot(2,3,6);imagesc(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1),[0 1]);
% title(sprintf('x_h (EPLL): PSNR %.1f SSIM %.3f',psnr(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))));
subplot(2,3,5);imagesc(X_hat(posM:posM+wM2-1, posN:wN2+posN-1),[0 1]);
title(sprintf('x_h (MMSE with one component): PSNR %.1f SSIM %.3f',psnr(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))));
subplot(2,3,6);imagesc(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1),[0 1]);
title(sprintf('x_h (MMSE with all components): PSNR %.1f SSIM %.3f',psnr(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))));
colormap gray
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 35 30]);
saveas(gcf,['fig/',img,'_parthd_',init],'fig'); % top-right
saveas(gcf,['png/',img,'_parthd_',init],'png');
tab_parthd=[psnr(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1));...
psnr(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1));...
psnr(xb(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(xb(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))];

figure()
posM=256;
posN=256;
posM2=posM/q;  
posN2=posN/q;
subplot(2,3,1);imagesc(x_orig(posM:posM+wM2-1, posN:wN2+posN-1),[0,1]);title('Partie de X');
subplot(2,3,2);imagesc(y_obs(posM2:posM2+wM-1, posN2:wN+posN2-1),[0 1]);
title('Partie de Y (connue)');
subplot(2,3,4);imagesc(xb(posM:posM+wM2-1, posN:wN2+posN-1),[0 1]);
title(sprintf('z: PSNR %.1f SSIM %.3f',psnr(xb(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(xb(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))));
% subplot(2,3,5);imagesc(X_hat(posM:posM+wM2-1, posN:wN2+posN-1),[0 1]);
% title(sprintf('x_h (MMSE): PSNR %.1f SSIM %.3f',psnr(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))));
% subplot(2,3,6);imagesc(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1),[0 1]);
% title(sprintf('x_h (EPLL): PSNR %.1f SSIM %.3f',psnr(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))));
subplot(2,3,5);imagesc(X_hat(posM:posM+wM2-1, posN:wN2+posN-1),[0 1]);
title(sprintf('x_h (MMSE with one component): PSNR %.1f SSIM %.3f',psnr(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))));
subplot(2,3,6);imagesc(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1),[0 1]);
title(sprintf('x_h (MMSE with all components): PSNR %.1f SSIM %.3f',psnr(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))));
colormap gray
tab_partbd=[psnr(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1));...
psnr(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1));...
psnr(xb(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(xb(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 35 30]);
saveas(gcf,['fig/',img,'_partbd_',init],'fig'); % bottom-right
saveas(gcf,['png/',img,'_partbd_',init],'png');

figure()
posM=256;
posN=1;
posM2=posM/q;  
posN2=1;
subplot(2,3,1);imagesc(x_orig(posM:posM+wM2-1, posN:wN2+posN-1),[0,1]);title('Partie de X');
subplot(2,3,2);imagesc(y_obs(posM2:posM2+wM-1, posN2:wN+posN2-1),[0 1]);
title('Partie de Y (connue)');
subplot(2,3,4);imagesc(xb(posM:posM+wM2-1, posN:wN2+posN-1),[0 1]);
title(sprintf('z: PSNR %.1f SSIM %.3f',psnr(xb(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(xb(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))));
% subplot(2,3,5);imagesc(X_hat(posM:posM+wM2-1, posN:wN2+posN-1),[0 1]);
% title(sprintf('x_h (MMSE): PSNR %.1f SSIM %.3f',psnr(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))));
% subplot(2,3,6);imagesc(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1),[0 1]);
% title(sprintf('x_h (EPLL): PSNR %.1f SSIM %.3f',psnr(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))));
subplot(2,3,5);imagesc(X_hat(posM:posM+wM2-1, posN:wN2+posN-1),[0 1]);
title(sprintf('x_h (MMSE with one component): PSNR %.1f SSIM %.3f',psnr(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))));
subplot(2,3,6);imagesc(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1),[0 1]);
title(sprintf('x_h (MMSE with all components): PSNR %.1f SSIM %.3f',psnr(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))));
colormap gray
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 35 30]);
saveas(gcf,['fig/',img,'_partbg_',init],'fig');  % bottom-left
saveas(gcf,['png/',img,'_partbg_',init],'png');
tab_partbg=[psnr(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1));...
psnr(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1));...
psnr(xb(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(xb(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))];

