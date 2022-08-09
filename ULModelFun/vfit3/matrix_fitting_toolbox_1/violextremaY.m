function [s_pass,g_pass,smin]=violextremaY(SERflag,wintervals,A,B,C,D,colinterch);

% [s_pass,g_pass]=violextremaY(SERflag,wintervals,A,B,C,D,colinterch);
%
% This file is part of the Matrix Fitting Toolbox, v1.
% Filename: violextremaY.m
% Package: Matrix_Fitting_Toolbox_1.zip.
% Programmed by B. Gustavsen. October 08, 2008.
%
%
% Purpose: Identify eigenvalue minima within given intervals   
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
%          g_pass: the smallest (most negative) eigenvalue encountered.   
%          smin  : the frequency where g_pass is located 

s_pass=[]; g_pass=[]; ss=[];

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
g_pass=1e16; smin=0;
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

  oldT0=[];
  for k=1:length(s_pass)
    Y = fitcalcABCDE(s_pass(k),diag(A),B,C,D,zeros(Nc));
    G=real(Y);
    if colinterch==0  
      EE(:,k)=eig(G);   
    else
      [T0,DD]=eig(G);
      T0=rot(T0); %minimizing phase angle of eigenvectors in least squares sense  
      [T0,DD]=intercheig(T0,oldT0,DD,Nc,k);  oldT0=T0;
      EE(:,k)=diag(DD);
    end  
  end

%====================================================
% Identifying violations, picking minima for s2:
%====================================================
s_pass_ind=zeros(1,length(s_pass));
for row=1:Nc
  if EE(row,1)<0
    s_pass_ind(1)=1;
  end
end  
%tell=0;
for k=2:length(s_pass)-1
  for row=1:Nc
    if EE(row,k)<0 %violation
      if (EE(row,k)<EE(row,k-1))&(EE(row,k)<EE(row,k+1))
        s_pass_ind(k)=1;
      end
    end
  end
end
s=[s s_pass(find(s_pass_ind==1))];
dum=min(EE);
[g_pass2,ind]=min(dum); smin2=s_pass(ind); %largest violation in interval
[g_pass,ind]=min([g_pass g_pass2]);
dums=[smin smin2];
smin=dums(ind);

g_pass=min(g_pass,min(min(EE)));

%[g_pass g_pass2 smin/(2*pi*i)]
end %for m=1:length(wintervals(:,1))

s_pass=s;


