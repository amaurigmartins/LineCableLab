function [s_pass,g_pass,smin]=violextremaS(SERflag,wintervals,A,B,C,D,colinterch);

% [s_pass,g_pass]=violextremaS(SERflag,wintervals,A,B,C,D,colinterch);
%
% This file is part of the Matrix Fitting Toolbox, v1.
% Filename: violextremaS.m
% Package: Matrix_Fitting_Toolbox_1.zip.
% Programmed by B. Gustavsen. October 08, 2008.
%
%
% Purpose: Identify singular value maxima within given intervals   
%
% Input:  sing_w:   array holding frequency where e-values of G cross zero line
%         intervals: array defining intervals of passivity violations by list of freq.pairs
%                       intervals(start_ind:end_ind,:).
%                       (start_ind_end_ind are indices into sing_w);  
%
%
% Output:  s_pass: array holding frequency samples where the minima occur.  
%                  For each interval, the minimum of all eigenvals are checked. There can
%                  at most be one minimum for each eigenvalue. Thus, each interval results in one or
%                  more entries in s2, with at most Nc entries. 
%          g_pass: the largest (most positive) singular value encountered. 

s_pass=[]; g_pass=[]; smin=[];

if length(wintervals)==0, s_pass=[]; return, end

if SERflag==1
  SERA=A; SERC=C; SERD=D; clear B;  
    
  Nc=length(SERD); N=length(SERA);
  A=zeros(Nc*N,1);
  for col=1:Nc
    A((col-1)*N+1:col*N)=SERA;
  end   
  B=zeros(Nc*N,Nc);
  for col=1:Nc
    B((col-1)*N+1:col*N,col)=ones(N,1);
  end
  C=zeros(Nc,Nc*N);
  for row=1:Nc
    for col=1:Nc 
      C(row,(col-1)*N+1:col*N)=squeeze(SERC(row,col,1:N)).';
    end    
  end
  D=SERD;
  
  A=diag(A);
end


s=[];
Nc=length(D);
tell=0;
for m=1:length(wintervals(:,1))

%For each interval, identify global minima for each e-value
  Nint=21; %n.o. internal frequency samples resolving each interval

  w1=wintervals(m,1);
  if wintervals(m,2)==1e16 %violation extends to infinite frequency 
    w2=2*pi*1e16; %tvilsom måte å spesifisere uendelig frekv. på
  else
    w2=wintervals(m,2);
  end
  s_pass1=i*linspace(w1,w2,Nint);
  if w1==0
    s_pass2=i*logspace(-8,log10(w2),Nint);
  else
    s_pass2=i*logspace(log10(w1),log10(w2),Nint);
  end
  s_pass=sort([s_pass1 s_pass2]); Nint=2*Nint;

  oldU=[];
  for k=1:length(s_pass)
    SS = fitcalcABCDE(s_pass(k),diag(A),B,C,D,zeros(Nc));
    if colinterch==0  
      EE(:,k)=svd(SS,0); 
    else 
      [U,S,V]=svd(SS,0);   
      [U,S,V]=interchsvd(U,oldU,S,V,Nc,k); oldU=U;   
      EE(:,k)=diag(S);   
    end ;
  end

%====================================================
% Identifying violations, picking minima for s2:
%====================================================
s_pass_ind=zeros(1,length(s_pass));
for row=1:Nc
  if EE(row,1)>1
    s_pass_ind(1)=1;
  end
end  
%tell=0;
for k=2:length(s_pass)-1
  for row=1:Nc
    if EE(row,k)>1 %violation
      if (EE(row,k)>EE(row,k-1))&(EE(row,k)>EE(row,k+1))
        s_pass_ind(k)=1;
        tell=tell+1; g_pass(tell)=EE(row,k); 
      end
    end
  end
end
s=[s s_pass(find(s_pass_ind==1))];


end %for m=1:length(wintervals(:,1))
s_pass=s;

g_pass=-1e16;
for k=1:length(s_pass);
  SS = fitcalcABCDE(s_pass(k),diag(A),B,C,D,zeros(Nc));
  dum(k)=max(svd(SS));
  %[g_pass ind]=max(g_pass,max(svd(SS)));
end
[g_pass ind]=max(dum);
smin=s_pass(ind);


