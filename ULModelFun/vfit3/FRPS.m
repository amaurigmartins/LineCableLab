function [SER,RPopts]=FRPS(SER,s,s2,s3,RPopts);

% [SER,RPopts]=FRPS(SER,s,s2,s3,RPopts);
%
% This file is part of the Matrix Fitting Toolbox, v1.
% Filename: FRPS.m
% Package: Matrix_Fitting_Toolbox_1.zip.
% Programmed by B. Gustavsen. October 08, 2008.
%
%      ======================================================
%      =   Routine: FRPS.m                                  = 
%      =   Version 1.0                                      =
%      =   Last revised: 08.10.2008                         = 
%      =   Bjorn Gustavsen                                  =
%      =   SINTEF Energy Research, N-7465 Trondheim, NORWAY =
%      =   bjorn.gustavsen@sintef.no                        =
%      ======================================================
%
%  PURPOSE : Modify elements in SER of rational approximation  to enforce passivity 
%            of S-parameter model at frequency samples in s2 and s3, such that the
%            perturbation of the model is minimized at samples in s. 
%
%  Reference: 
%     B. Gustavsen and A. Semlyen, "Fast passivity assessment for S-parameter rational
%     models via a half-size test matrix", 
%     IEEE Trans. Microwave Theory and Techniques, accepted for publication. 

auxflag=RPopts.auxflag;
weightfactor=RPopts.weightfactor;
weightparam=RPopts.weightparam;
TOL=RPopts.TOLGD;


SERA=SER.poles;
SERAA=SER.A;
SERBB=SER.B;
SERCC=SER.C;
[m,n]=size(SERA);
if m<n
  SERA=SERA.'; %Må ha kolonne-vektor
end  
SERC=SER.R; 
SERD=SER.D;
SERE=SER.E;


if ~isfield(RPopts,'H')
  RPopts.H=[]; 
  RPopts.oldDflag=-1;
end

%d=eig(SERD);eigD=d;
[U,S,V]=svd(SERD,0);
if any(diag(S)>1)
  Dflag=1; %Will perturb D-matrix
  [VD,eigD]=eig(SERD); invVD=VD^(-1); eigD=diag(eigD);
else
  Dflag=0;  
end

SERCnew=SERC; SERDnew=SERD; SEREnew=SERE;

N=length(SERA);

bigB=[];
bigc=[];
Ns=length(s);
Ns2=length(s2);
Nc=length(SERD);
%!!Nc2=Nc*Nc;
Nc2=Nc*(Nc+1)/2;
I=eye(Nc);
M2mat=[];


%=======================================================
% Finding out which poles are complex :
%=======================================================
% LS problem:
cindex=zeros(1,N);
for m=1:N 
  if imag(SERA(m))~=0  
    if m==1 
      cindex(m)=1;
    else
      if cindex(m-1)==0 | cindex(m-1)==2
        cindex(m)=1; cindex(m+1)=2; 
      else
        cindex(m)=2;
      end
    end 
  end
end



    
if isempty(RPopts.H) %ÆÆÆÆÆÆÆÆÆÆÆÆÆÆÆÆÆÆÆÆÆÆÆ

if RPopts.outputlevel==1, ,disp('    Building system equation (once)...'), end 
%========================================
% LOOP FOR LEAST SQUARES PROBLEM:
%========================================
if (Dflag)==1    
  bigArow=zeros(Nc2,Nc*(N+1));
  H=zeros(Nc*(N+1));
  Mmat=zeros(Nc2,Nc*(N+1));
else    
  bigArow=zeros(Ns*Nc2,Nc*(N));
  H=zeros(Nc*(N));  
  Mmat=zeros(Nc2,Nc*(N));  
end  
for m=1:N
  R=(SERC(:,:,m));
  if cindex(m)==0      %real pole  
    R=R;
  elseif cindex(m)==1  %complex pole, 1st part
    R=real(R);
  else  
    R=imag(R);
  end  
  [V,D]=eig(R);
  bigV(1:Nc,(m-1)*Nc+1:m*Nc)=V;
  biginvV(1:Nc,(m-1)*Nc+1:m*Nc)=V^(-1);  
  bigD(:,m)=diag(D);
end  

%Beregning av kolonne skalering:
if Dflag==1,Escale=zeros(N+1,1);else,Escale=zeros(N,1);end
for m=1:N
  if cindex(m)==0      %real pole
    Escale(m)=1/norm( 1./(s-SERA(m)) );
  elseif cindex(m)==1  %complex pole, 1st part
    Escale(m)  =1/norm( 1./(s-SERA(m)) + 1./(s-SERA(m)') );
  else  
    Escale(m)  =1/norm( i./(s-SERA(m)') - i./(s-SERA(m)) ); 
  end
end      
if Dflag==1,Escale(N+1)=1./sqrt(Ns); end

%Escale=1+0.*Escale;
if Dflag==1, dum=N+1;else,dum=N;end
Escale2=zeros(Nc*dum,1);
for m=1:dum
  Escale2((m-1)*Nc+1:Nc*dum)=Escale(m);
end 
RPopts.Escale=Escale;
RPopts.Escale2=Escale2;
%Creatiung indices to lower triangle:
indices=zeros(Nc2,1);
tell=0;
for col=1:Nc
  for row=col:Nc
    tell=tell+1;  
    indices(tell)=(col-1)*Nc+row;
  end  
end
for k=1:Ns
    %k
      sk=s(k);
      % Calculating matrix Mmat (coefficients for LS-problem)
      tell=0;
      offs=0;
      %Yfit=Yfitcalc(sk,SERA,SERC,SERD,SERE);   
      Yfit=fitcalcABCDE(sk,diag(SERAA),SERBB,SERCC,SERD,SERE);   
      %y=SERD(row,col)+sk*SERE(row,col);
      %y=y+sum(squeeze(SERC(row,col,1:N))./(s(k)-SERA(1:N)));
      if weightparam==1
        weight=ones(Nc,Nc);
      elseif weightparam==2
        weight=1./abs(Yfit);          
      elseif weightparam==3
        weight=1./sqrt(abs(Yfit));
      elseif weightparam==4 
        weight=ones(Nc,Nc)./norm(abs(Yfit));   
      elseif weightparam==5 
        weight=ones(Nc,Nc)./sqrt(norm(abs(Yfit)));          
      else
        'ERROR!!!!!!!!!'
        weight=1;
      end
  
  for m=1:N
    V=(bigV(:,(m-1)*Nc+1:m*Nc));
    invV=V^(-1);%squeeze(biginvV(:,(m-1)*Nc+1:m*Nc)); 
    if cindex(m)==0 %real pole
      dum=Escale(m)/(sk-SERA(m));   
    elseif cindex(m)==1  %complex pole, 1st part
      dum=Escale(m)*(1/(sk-SERA(m)) + 1/(sk-SERA(m)') );   
    else
      dum=Escale(m)*(i/(sk-SERA(m)') - i/(sk-SERA(m)) ); 
    end 
    
    for egenverdi=1:Nc
      tell=0;
      gamm=V(:,egenverdi)*invV(egenverdi,:);
      gammx=gamm.*weight.*dum;
      gammy=zeros(Nc2,1);
      Mmat(:,offs+egenverdi)=gammx(indices);
     end %for egenverdi=1:Nc
    offs=offs + (Nc);
  end %for m=1:N
 
  if Dflag==1
    for egenverdi=1:Nc %Egenverdi til residue-matrise
      gamm=Escale(N+1)*VD(:,egenverdi)*invVD(egenverdi,:); %Ytre-produkt 
      tell=0;
      for row=1:Nc
        for col=row:Nc
          tell=tell+1;  
          if row==col
            faktor=weight(row,col);  
          else
            faktor=2*weight(row,col);                
          end 
          Mmat(tell,offs+egenverdi)=gamm(row,col)*faktor; 
        end
      end  
    end  
  end  
  
  XX=[real(Mmat);imag(Mmat)]; 
  H=H+XX.'*XX;
  
end %for k=1:Ns

%========================================
% INTRODUCING SAMPLES OUTSIDE LS REGION: ONE SAMPLE PER POLE (s4)
%========================================
if auxflag==1
s4=[];    
tell=0;
for m=1:length(SERA)
  if cindex(m)==0 %real pole
    if ( ( abs(SERA(m))>s(Ns)/j ) | ( abs(SERA(m))<s(1)/j ) )
      tell=tell+1;
      s4(tell)=j*abs(SERA(m));
    end
  elseif cindex(m)==1 %complex pole, first part
    if ( ( abs(imag(SERA(m)))>s(Ns)/j ) | ( abs(imag(SERA(m)))<s(1)/j ) )
      tell=tell+1;
      s4(tell)=j*abs(imag(SERA(m)));
    end
  end
end
Ns4=length(s4);

bigA2=zeros(Ns4*Nc2,Nc*(N+Dflag));
 
for k=1:1:Ns4
      sk=s4(k);
      % Calculating matrix Mmat (coefficients for LS-problem)
      tell=0;
      offs=0;
      Yfit=fitcalcPRE(sk,SERA,SERC,SERD,SERE);
      %y=SERD(row,col)+sk*SERE(row,col);
      %y=y+sum(squeeze(SERC(row,col,1:N))./(s(k)-SERA(1:N)));
      if weightparam==1
        weight=ones(Nc,Nc);
      elseif weightparam==2
        weight=1./abs(Yfit);  
      elseif weightparam==3
        weight=1./sqrt(abs(Yfit)); 
      elseif weightparam==4
        weight=ones(Nc,Nc)./norm(abs(Yfit));  
      elseif weightparam==5 
        weight=ones(Nc,Nc)./sqrt(norm(abs(Yfit)));     
      else
        'ERROR!!!!!!!!!'
        weight=1;
      end
      weight=weight*weightfactor;
      
  for m=1:N
    V=squeeze(bigV(:,(m-1)*Nc+1:m*Nc));
    invV=V^(-1);%squeeze(biginvV(:,(m-1)*Nc+1:m*Nc)); 
    if cindex(m)==0 %real pole
      dum=Escale(m)/(sk-SERA(m));   
    elseif cindex(m)==1  %complex pole, 1st part
      dum=Escale(m)*(1/(sk-SERA(m)) + 1/(sk-SERA(m)') );   
    else
      dum=Escale(m)*(i/(sk-SERA(m)') - i/(sk-SERA(m)) ); 
    end 
    
    for egenverdi=1:Nc
      tell=0;
      gamm=V(:,egenverdi)*invV(egenverdi,:);     
      for row=1:Nc
        for col=row:Nc
          tell=tell+1;  
          if row==col
            faktor=weight(row,col);  
          else
            faktor=2*weight(row,col);                
          end   
          Mmat(tell,offs+egenverdi)=gamm(row,col)*faktor*dum;    
        end  
      end  
    end %for egenverdi=1:Nc
    offs=offs + (Nc);
  end %for m=1:N
 
  if Dflag==1
    for egenverdi=1:Nc %Egenverdi til residue-matrise
      gamm=Escale(N+1)*VD(:,egenverdi)*invVD(egenverdi,:); %Ytre-produkt 
      tell=0;
      for row=1:Nc
        for col=row:Nc
          tell=tell+1;
          if row==col
            faktor=weight(row,col);  
          else
            faktor=2*weight(row,col);                
          end    
          Mmat(tell,offs+egenverdi)=gamm(row,col)*faktor; 
        end
      end  
    end  
  end  
  
  H=H+ real(Mmat.')*real(Mmat) + imag(Mmat.')*imag(Mmat);
end %for k=1:Ns


end %if auxflag==1


RPopts.H=H; clear H
RPopts.bigV=bigV;
RPopts.biginvV=biginvV;
if RPopts.outputlevel==1, disp('    Done'), end 

else %if isempty(RPopts)
  bigV=RPopts.bigV;  
  biginvV=RPopts.biginvV;
  if Dflag~=RPopts.oldDflag 
    RPopts.H=RPopts.H(1:Nc*(N+Dflag),1:Nc*(N+Dflag));
    RPopts.Escale=RPopts.Escale((N+Dflag));   
  end  

end %if isempty(RPopts) 

Escale=RPopts.Escale;
Escale2=RPopts.Escale2;

Mmat2=zeros(Nc,Nc*(N+Dflag));

%========================================
% LOOP FOR CONSTRAINT PROBLEM, TYPE #1 (violating eigenvalues in s2):
%========================================

viol_S=[];
viol_D=[];

for k=1:Ns2
  sk=s2(k);

%   for row=1:Nc
%     for col=1:Nc
%       Y(row,col)=SERD(row,col)+sk*SERE(row,col);
%       Y(row,col)=Y(row,col)+sum(squeeze(SERC(row,col,1:N))./(sk-SERA(1:N)));
%     end
%   end 
  Y=fitcalcABCDE(sk,diag(SERAA),SERBB,SERCC,SERD,SERE);

    [U,S,V]=svd(Y,0);
    SS=diag(S);
    Qa=0.5*V^(-1);
    Qb=0.5*U^(-1);   
    Ta=V;
    Tb=U; 

  if max(diag(S))>1  %any violations
    % Calculating matrix M2mat; matrix of partial derivatives :
    tell=0;
    offs=0;
  
  for m=1:N    
    if cindex(m)==0 %real pole   
      VV=bigV(:,(m-1)*Nc+1:m*Nc);
      for egenverdi=1:Nc 
       tell=0;
       p=VV(:,egenverdi);
       for row=1:Nc %singular values
          tell=tell+1;
          k1=(Qb(row,:)*p)*(p.'*Ta(:,row)); 
          Mmat2(tell,offs+egenverdi)=2*k1/(sk-SERA(m));
        end   
      end  
      offs=offs + (Nc);  
    elseif cindex(m)==1  %complex pole, 1st part       
      VV1=      bigV(:,(m-1)*Nc+1:m*Nc); 
      VV2=      bigV(:,(m+1-1)*Nc+1:(m+1)*Nc); 
      f=1/(sk-SERA(m)) + 1/(sk-SERA(m)');  
      g=i/(sk-SERA(m)) - i/(sk-SERA(m)');
      for egenverdi=1:Nc 
        tell=0;
        p1=VV1(:,egenverdi);
        p2=VV2(:,egenverdi);     
        for row=1:Nc %singular values
          tell=tell+1;
          k1=2*f*(Qb(row,:)*p1)*(p1.'*Ta(:,row)); 
          k2=2*g*(Qb(row,:)*p2)*(p2.'*Ta(:,row));   
          Mmat2(tell,offs+egenverdi)=k1; 
          Mmat2(tell,offs+egenverdi+Nc)=k2;     
        end  
      end %for egenverdi=1:Nc 
      offs=offs + 2*(Nc);
    end 
  end %for m=1:N
  
   if Dflag==1       
    for egenverdi=1:Nc 
      tell=0;
      gamm=VD(:,egenverdi)*invVD(egenverdi,:);     
      for row=1:Nc %singular values
         k1= Qb(row,:)*gamm*Ta(:,row);  
         tell=tell+1;
         Mmat2(tell,offs+egenverdi)=2*k1;
      end  
    end %for egenverdi=1:Nc
   end  
    B=Mmat2;
    dels=SS;   
    for n=1:Nc   %ustabilitet?
      if dels(n)>1
        bigB=[bigB;B(n,:)];
        bigc=[bigc;-TOL+1-dels(n)];     
        viol_S=[viol_S;dels(n)];
      end
    end
          
  end %if max(real(Z))>0  %any violations

end %for k=1:Ns2

%bigB=real(bigB);

%========================================
% LOOP FOR CONSTRAINT PROBLEM, TYPE #2: (all eigenvalues in s3):
%========================================
Ns3=length(s3);
for k=1:Ns3
  sk=s3(k);

  Y=fitcalcABCDE(sk,diag(SERAA),SERBB,SERCC,SERD,SERE);

    [U,S,V]=svd(Y,0);
    SS=diag(S);
    Qa=0.5*V^(-1);
    Qb=0.5*U^(-1);   
    Ta=V;
    Tb=U; 

  if max(diag(S))>1  %any violations
    % Calculating matrix M2mat; matrix of partial derivatives :
    tell=0;
    offs=0;
  
  for m=1:N    
    if cindex(m)==0 %real pole   
      VV=bigV(:,(m-1)*Nc+1:m*Nc);
      for egenverdi=1:Nc 
       tell=0;
       p=VV(:,egenverdi);
       for row=1:Nc %singular values
          tell=tell+1;
          k1=(Qb(row,:)*p)*(p.'*Ta(:,row)); 
          Mmat2(tell,offs+egenverdi)=2*k1/(sk-SERA(m));
        end   
      end  
      offs=offs + (Nc);  
    elseif cindex(m)==1  %complex pole, 1st part       
      VV1=      bigV(:,(m-1)*Nc+1:m*Nc); 
      VV2=      bigV(:,(m+1-1)*Nc+1:(m+1)*Nc); 
      f=(1/(sk-SERA(m)) + 1/(sk-SERA(m)'));  
      g=(i/(sk-SERA(m)) - i/(sk-SERA(m)'));
      for egenverdi=1:Nc 
        tell=0;
        p1=VV1(:,egenverdi);
        p2=VV2(:,egenverdi);     
        for row=1:Nc %singular values
          tell=tell+1;
          k1=2*f*(Qb(row,:)*p1)*(p1.'*Ta(:,row)); 
          k2=2*g*(Qb(row,:)*p2)*(p2.'*Ta(:,row));   
          Mmat2(tell,offs+egenverdi)=k1; 
          Mmat2(tell,offs+egenverdi+Nc)=k2;     
        end  
      end %for egenverdi=1:Nc 
      offs=offs + 2*(Nc);
    end 
  end %for m=1:N
  
   if Dflag==1       
    for egenverdi=1:Nc 
      tell=0;
      gamm=Escale(m+1)*VD(:,egenverdi)*invVD(egenverdi,:);     
      for row=1:Nc %singular values
         k1= Qb(row,:)*gamm*Ta(:,row);  
         tell=tell+1;
         Mmat2(tell,offs+egenverdi)=2*k1;
      end  
    end %for egenverdi=1:Nc
   end  
    B=Mmat2;
    dels=SS;   
    for n=1:Nc   %ustabilitet?
      %if dels(n)>1
        bigB=[bigB;B(n,:)];
        bigc=[bigc;-TOL+1-dels(n)];    
        viol_S=[viol_S;dels(n)];      
      %end
    end
          
  end %if max(real(Z))>0  %any violations

end %for k=1:Ns3

%Legger til constraint for evt. sing.-verdier i D >1
 if Dflag==1   
  %!![U,S,V]=svd(Escale(m+1)*SER.D,0);   
  [U,S,V]=svd(SER.D,0);  
  SS=diag(S);
  Qa=0.5*V^(-1);
  Qb=0.5*U^(-1);   
  Ta=V;
  Tb=U;
  dumB=zeros(sum(SS>1),Nc*N+Nc);
  dumc=zeros(sum(SS>1),1);
  for egenverdi=1:Nc 
    tell=0;
    gamm=VD(:,egenverdi)*invVD(egenverdi,:);     
    for row=1:Nc %singular values
      if SS(row)>1  
        koko=[Qa(row,:) Qb(row,:)]*[Ta(:,row);Tb(:,row)]; 
        k1= Qb(row,:)*gamm*Ta(:,row); 
        k2= Qa(row,:)*gamm*Tb(:,row);     
        tell=tell+1;
        dumB(tell,Nc*N+egenverdi)=(1/koko)*( k1 + k2 );  
        dumc(tell)=-TOL+1-SS(row);
        viol_S=[viol_S; SS(row)];     
        viol_D=[viol_D; SS(row)];  
      end  
      %end  
    end  
  end %for egenverdi=1:Nc
   bigB=[bigB;dumB];
   bigc=[bigc;dumc];  
 end


if length(bigB)==0 
  return   %No passivity violations
end


c=bigc;


bigB=real(bigB);
for col=1:length(RPopts.H)
  if length(bigB)>0
    bigB(:,col)=bigB(:,col).*RPopts.Escale2(col);
  end
end  

%bigB=real(bigB)


ff=zeros(length(RPopts.H),1);

clear bigA;
%disp('Starting quadprog...')

warning off


if RPopts.solver=='QUADPROG'
  [dx,lambda]=quadprog(RPopts.H,ff,bigB,bigc);  
elseif RPopts.solver=='CPLEX   '  

  c=0*ff;
  H0=c.';

  %F=sparse(H);
  A=-bigB;

  x_0=zeros(length(H(:,1)));
  b_U=bigc;
  b_L = -inf*ones(length(b_U),1);
  x_L = -inf*ones(length(x_0),1);
  x_U = inf*ones(length(x_0),1);

  Prob = qpAssign(F,c,A,b_L,b_U,x_L,x_U,x_0,'dust');
  Prob.MIP.cpxControl.QPMETHOD = 4;
  Prob.MIP.cpxControl.BARALG = 2;
  PriLev = 0;
  Prob.PriLevOpt = 0;
  [x, slack, v, rc, f_k] = cplex(c, A, x_L, x_U, b_L, b_U, ...
           Prob.MIP.cpxControl, [], PriLev, Prob, [], [], [], [], ...
           [], [], F); 
  dx=x;

end
%disp(['Elapsed time: ' num2str(toc) ' sec'])
warning on

dx=dx.*RPopts.Escale2(1:Nc*(N+Dflag));


%Oppdaterer egenverdiene og SERC:
for m=1:N
  if cindex(m)==0 %real pole
    D1=diag(dx((m-1)*Nc+1:m*Nc));  
    SERCnew(:,:,m)=SERCnew(:,:,m) +( bigV(:,(m-1)*Nc+1:m*Nc)*D1*biginvV(:,(m-1)*Nc+1:m*Nc) );  
  elseif cindex(m)==1  %complex pole, 1st part
    GAMM1=bigV(:,(m-1)*Nc+1:m*Nc);
    GAMM2=bigV(:,(m+1-1)*Nc+1:(m+1)*Nc);
    invGAMM1=biginvV(:,(m-1)*Nc+1:m*Nc);
    invGAMM2=biginvV(:,(m+1-1)*Nc+1:(m+1)*Nc);     
      
    D1=diag(dx((m-1)*Nc+1:m*Nc));          
    D2=diag(dx((m+1-1)*Nc+1:(m+1)*Nc));  
    R1=real(SERC(:,:,m));
    R2=imag(SERC(:,:,m));    
    R1new=R1 +GAMM1*D1*invGAMM1; 
    R2new=R2 +GAMM2*D2*invGAMM2; 
    SERCnew(:,:,m)  =R1new+j*R2new;
    SERCnew(:,:,m+1)=R1new-j*R2new;
  end  
end

if Dflag==1
  DD=diag(dx(N*Nc+1:(N+1)*Nc));
  SERDnew=SERDnew+VD*DD*invVD;
end

%==========================
SERDnew=(SERDnew+SERDnew.')/2;

for m=1:N
  SERCnew(:,:,m)=(SERCnew(:,:,m)+SERCnew(:,:,m).')/2;
end  
%==========================
SER.R=SERCnew;
SER.D=SERDnew;
SER.E=zeros(Nc);
[SER]=pr2ss(SER);  

RPopts.oldDflag=Dflag;




end %function FRPS 


