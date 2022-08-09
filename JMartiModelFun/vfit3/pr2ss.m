function [SER]=pr2ss(SER);

% [SER]=pr2ss(SER);
%
% This file is part of the Matrix Fitting Toolbox, v1.
% Filename: pr2ss.m
% Package: Matrix_Fitting_Toolbox_1.zip.
% Programmed by B. Gustavsen. October 08, 2008.

R=SER.R;
poles=SER.poles;
Nc=length(SER.D);
N=length(R(1,1,:));
C=zeros(Nc,Nc*N);
A=sparse(Nc*N,1);
B=sparse(Nc*N,Nc);
poles=diag(diag(poles)); %Blir kolonne-vektor
for m=1:N   
  Rdum=squeeze(R(:,:,m));  
  for n=1:Nc
    ind=(n-1)*N+m;
    C(:,ind)=Rdum(:,n);
  end
end
for n=1:Nc
  A((n-1)*N+1:n*N)=poles;
  B((n-1)*N+1:n*N,n)=ones(N,1);  
end
A=diag(A);

SER.A=A;
SER.B=B;
SER.C=C;

end