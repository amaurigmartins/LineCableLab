function [SER1]=reducecmplx(SER); 
%
% [SER1]=reducecmplx(SER); 
% 
% Processing of model before time domain simulation, for improved speed.
% The Model size is reduced by throwing out the second part of all complex pairs
% and multiplying rows in B associated with first part by factor 2.  
%
%  Download site:
%  http://www.energy.sintef.no/Produkt/VECTFIT/index.asp 
%
%  17.03.2013. Bjorn Gustavsen, SINTEF Energy Research, Norway.


N=length(SER.A);
Nc=length(SER.D);


  
A=SER.A; B=SER.B; C=SER.C; D=SER.D; E=SER.E; 

%##################################
N=length(SER.A);
cindex=zeros(1,N);
for m=1:N 
  if imag(SER.A(m,m))~=0  
    if m==1 
      cindex(m)=1;
    else
      if cindex(m-1)==0 || cindex(m-1)==2
        cindex(m)=1; cindex(m+1)=2; 
      else
        cindex(m)=2;
      end
    end 
  end
end
ind=find(cindex==2);
C(:,ind-1)=C(:,ind-1).*2; 
C(:,ind)=[];A(ind,:)=[];A(:,ind)=[];B(ind,:)=[];
N=length(A); 

SER1.A=A;
SER1.B=B;
SER1.C=C;
SER1.D=D;
SER1.E=E;


end  %function 

