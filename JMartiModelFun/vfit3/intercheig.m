function [V,D]=intercheig(V,oldV,D,Nc,fstep);   

% [V,D]=intercheig(V,oldV,D,Nc,fstep);
%
% This file is part of the Matrix Fitting Toolbox, v1.
% Filename: intercheig.m
% Package: Matrix_Fitting_Toolbox_1.zip.
% Programmed by B. Gustavsen. October 08, 2008.
%
% Interchanging e-vectors and corresponding e-values so that they become smooth
% functions of frequency. Considers the old eigenvectors one-by-one and finds out by which of 
% the new ones the largest dot product is attained. 
% [V,D]=eig(A)


if fstep>1 
    
  UGH=abs(real(oldV'*V));  
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

     hjelp=V(:,ii);
     V(:,ii)=V(:,ilargest);
     V(:,ilargest)=hjelp;

     hjelp=D(ii,ii);
     D(ii,ii)=D(ilargest,ilargest);
     D(ilargest,ilargest)=hjelp;
     
     dum=UGH(:,ii);
     UGH(:,ii)=UGH(:,ilargest); UGH(:,ilargest)=dum;
  end 

%C Finding out whether the direction of e-vectors are 180deg. in error. Is done by comparing
%C sign of dotproducts of e-vectors, for new and old V-matrix :
  for ii=1:Nc
    dotprod=oldV(:,ii)'*V(:,j);        
    if sign(real(dotprod))<0 
      V(:,ii)=-V(:,ii);
    end
  end

end %
