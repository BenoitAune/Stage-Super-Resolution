function [Yd, we]=MAP_niknejad(Y,Ya,sig,gamma,H,Subt,Sub)

%input: Y contains the patches of the reconstructed image after the knn 
%       Ya contains the patches of the observed image after the knn 
%       sig is the variance of noise
%       gamma is a weights coeficients for the image reconstruction
%       H is the blur operator
%       Sub & Subt are the sub sampling operators. (Sub to shrink & Subt to extend)
%output: Yd contains the patches of the new reconstructed image
%       we are the weights used to reconstruct the image

    Yd=zeros(size(Y));
    epsilon=0.1;
    mu=mean(Y,2);
    for iter=1:1
        dimy=size(Y,2);
        Sc=zeros(size(Y,1),size(Y,1));
    for i=1:dimy
        y=Y(:,i);
        y=y-mu;
        covy=(y)*(y)';
        Sc=covy+Sc;   
    end
    covY=Sc./(dimy-1);
    covY=covY+epsilon*eye(size(covY));  % covariance computation
    [U ,S]=eig(covY);                   % to compute easily the inverse of the covariance
    
    Si=zeros(size(S));
    m=size(S,1);
    for i=1:m
            S(i,i)=S(i,i)-0;
            Si(i,i)=1/(S(i,i)+epsilon);
    end
    S=S+epsilon*eye(size(S));
    dets=1;
    for i=1:m
        dets=dets*S(i,i);
    end
    %U*Si*U' is the inverse
    las=10; % is a factor to increase the denoising 
    % we use the MAP
     HtH=H*Subt*Sub*H;
     gauche =(HtH+las*sig^2*(U*Si*U'));
    for i=1:size(Y,2)
      Yr=Subt*Ya(:,i);
      Hty=H*Yr;
      droite=Hty+(las*sig^2*(U*Si*U')*mu);      
      Yd(:,i)=gauche\droite;
       l=(las*sig^2*(Yd(:,i)-mu)'*(U*Si*U')*(Yd(:,i)-mu)+las*sig^2*log(dets))*(gamma);
       we(i,1)=exp(-l);
    end   
    end
    
    end