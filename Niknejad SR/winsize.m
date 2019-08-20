function [sr,sc,er,ec]=winsize(rowp,colp,N,dim,dim2)
difr=0;
difc=0;
if rowp-N/2>0
    sr=rowp-N/2;
else
    difr=N/2-rowp;
    sr=1;
end
if ~(rowp+N/2>dim)
    er=rowp+N/2+difr;
else
    er=dim;
    difr=rowp+N/2-dim;
    sr=sr-difr;
end
if colp-N/2>0
    sc=colp-N/2;
else
    sc=1;
    difc=N/2-colp;
end
if ~(colp+N/2>dim2)
    ec=colp+N/2+difc;
else
    ec=dim2;
    difc=colp+N/2-dim2;
    sc=sc-difc;
end
    
    
    
    
    