function [Yd, we]=lowrcovsimp2(Y,Ya,sig,gamma,G)

  Yd=zeros(size(Y));
    las=1;
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
    covY=covY+epsilon*eye(size(covY));
    [U ,S]=eig(covY);
    
    Si=zeros(size(S));
    m=size(S,1);
    for i=1:m
            S(i,i)=S(i,i)-0;
            Si(i,i)=1/(S(i,i)+epsilon);
    end
    S=S+epsilon*eye(size(S));
    
      
    Ym=zeros(size(Y));
    for i=1:size(Y,2)
        Ym(:,i)=mu;
    end
    X=U'*Ya+1*(las*(sig^2)*Si)*U'*Ym;
    for i=1:size(S,1)
    X(i,:)=(1/(1+las*(sig^2)*Si(i,i)))*X(i,:);
    end  
    Yd=U*X;
    
    
    icovY=inv(covY);
    dets=1;
    for i=1:m
        dets=dets*S(i,i);
    end
    we=ones(size(Y,2),1);

    for i=1:size(Y,2)
      
      wwi=inv(G*G+sig^2*icovY);
      Yd(:,i)=wwi*(G*Ya(:,i)+sig^2*icovY*mu);
      l=((Yd(:,i)-Ya(:,i))'*(Yd(:,i)-Ya(:,i))+las*sig^2*(Yd(:,i)-mu)'*(U*Si*U')*(Yd(:,i)-mu)+las*sig^2*log(dets))*gamma;
      we(i,1)=exp(-l);
      
    end
    
    end