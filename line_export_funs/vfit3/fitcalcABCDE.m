function Yfit=fitcalcABCDE(sk,A,B,C,D,E);

% Yfit=fitcalcABCDE(sk,A,B,C,D,E);
%
% This file is part of the Matrix Fitting Toolbox, v1.
% Filename: fitcalcABCDE.m
% Package: Matrix_Fitting_Toolbox_1.zip.
% Programmed by B. Gustavsen. October 08, 2008.

Nc=length(D);
N=length(A);
I=ones(N,1);
% for m=1:N
%   C(:,m)=C(:,m)./(sk-A(m));
% end

% for n=1:Nc
%   C(n,:)=C(n,:)./(sk-A.');
% end  

dum=repmat(1./(sk-A),1,Nc);
C=C.*dum.';

Yfit=C*B +D +sk*E;

