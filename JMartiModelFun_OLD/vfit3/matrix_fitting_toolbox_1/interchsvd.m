function [U,S,V]=interchsvd(U,oldU,S,V,Nc,fstep);  

% [U,S,V]=interchsvd(U,oldU,S,V,Nc,fstep);  
%
% This file is part of the Matrix Fitting Toolbox, v1.
% Filename: interchsvd.m
% Package: Matrix_Fitting_Toolbox_1.zip.
% Programmed by B. Gustavsen. October 08, 2008.
%
% Interchanging singular values and -vectors so that they become smooth
% functions of frequency. Considers the old singular vectors one-by-one and finds out by which of 
% the new ones the largest dot product is attained. 
% [U,S,V]=svd(A)


% Translation: , hjelp=help,

if fstep>1 

  UGH=abs(real(oldU'*U));   
  for ii=1:Nc 
    ilargest=0;
    rlargest=0;
    for j=1:Nc
      dotprod=UGH(ii,j);
      if dotprod>rlargest           
        rlargest=abs(real(dotprod));
        ilargest=j;       
      end
    end
    dot(ii)=rlargest;
    ind(ii)=ii;       
    taken(ii)=0;                                       
  end %

%C Sorting inner products abs(realdel) in descending order :
  for ii=1:Nc  
    for j=1:Nc-1
      if(dot(j)<dot(j+1));
        hjelp(1)=dot(j+1);
        ihjelp  =ind(j+1);
        dot(j+1)=dot(j);
        ind(j+1)=ind(j);
        dot(j)=hjelp(1);
        ind(j)=ihjelp;
      end
    end
  end

%C Doing the interchange in a prioritized sequence :
  for l=1:Nc %
    ii=ind(l);
    ilargest=0;
    rlargest=0;

    for j=1:Nc 
      if taken(j)==0   
        dotprod=UGH(ii,j);
        if dotprod>rlargest          
          rlargest=abs(real(dotprod));
          ilargest=j;       
        end
      end 
    end 

    taken(ii)=1;

     hjelp=U(:,ii);
     U(:,ii)=U(:,ilargest);
     U(:,ilargest)=hjelp;

     hjelp=S(ii,ii);
     S(ii,ii)=S(ilargest,ilargest);
     S(ilargest,ilargest)=hjelp;
     
     hjelp=V(:,ii);
     V(:,ii)=V(:,ilargest);
     V(:,ilargest)=hjelp;
     dum=UGH(:,ii);
     UGH(:,ii)=UGH(:,ilargest); UGH(:,ilargest)=dum;
  end 

%C Finding out whether the direction of sing-vectors are 180deg. in error. Is done by comparing
%C sign of dotproducts of e-vectors, for new and old U-matrix :
  for ii=1:Nc

    dotprod=(oldU(:,ii)'*U(:,ii));      
    if sign(real(dotprod))<0 
      U(:,ii)=-U(:,ii);
      V(:,ii)=-V(:,ii);
    end    
  end

end %

