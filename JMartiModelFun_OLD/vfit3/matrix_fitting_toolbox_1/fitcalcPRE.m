function Yfit=fitcalcPRE(s,SERA,SERC,SERD,SERE);

% Yfit=fitcalcPRE(s,SERA,SERC,SERD,SERE);
%
% This file is part of the Matrix Fitting Toolbox, v1.
% Filename: fitcalcPRE.m
% Package: Matrix_Fitting_Toolbox_1.zip.
% Programmed by B. Gustavsen. October 08, 2008.

Ns=length(s);
Nc=length(SERD);
N=length(SERA);
Yfit=zeros(Nc,Nc,Ns);
for k=1:Ns
  tell=0;
  for row=1:Nc
    for col=1:Nc
      tell=tell+1;
      Y(row,col)=SERD(row,col)+s(k)*SERE(row,col);
      Y(row,col)=Y(row,col)+sum(squeeze(SERC(row,col,1:N))./(s(k)-SERA(1:N)));   
    end
  end
  Yfit(1:Nc,1:Nc,k)=Y;
end

end 

