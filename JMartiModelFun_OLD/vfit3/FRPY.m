function [SER,RPopts]=FRPY(SER,s,s2,s3,RPopts);

% [SER,RPopts]=FRPY(SER,s,s2,s3,RPopts);
%
% This file is part of the Matrix Fitting Toolbox, v1.
% Filename: FRPY.m
% Package: Matrix_Fitting_Toolbox_1.zip.
% Programmed by B. Gustavsen. October 08, 2008.
%
%      ======================================================
%      =   Routine: FRPY.m                                  = 
%      =   Version 1.0                                      =
%      =   Last revised: 27.09.2009                         = 
%      =   Bjorn Gustavsen                                  =
%      =   SINTEF Energy Research, N-7465 Trondheim, NORWAY =
%      =   bjorn.gustavsen@sintef.no                        =
%      ======================================================
%
%  PURPOSE : Modify elements in SER of rational approximation  to enforce passivity 
%            of Y-parameter model at frequency samples in s2 and s3, such that the
%            perturbation of the model is minimized at samples in s. 
%
%  Reference: 
%     B. Gustavsen, "Fast passivity enforcement for S-parameter models by perturbation
%     of residue matrix eigenvalues", 
%     IEEE Trans. Advanced Packaging, accepted for publication. 

auxflag=RPopts.auxflag;
weightfactor=RPopts.weightfactor;
weightparam=RPopts.weightparam;
if isfield(RPopts,'weight'),
  bigweight=RPopts.weight;
end
TOLE=RPopts.TOLE;
TOL=RPopts.TOLGD;

SERA=SER.poles;
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
  RPopts.oldEflag=-1;  
end

d=eig(SERD);eigD=d;
if any(d<0)
  Dflag=1; %Will perturb D-matrix
  [VD,eigD]=eig(SERD); invVD=VD^(-1); eigD=diag(eigD);
else
  Dflag=0;  
end

e=eig(SERE);eigE=e;
if any(e<0)
  Eflag=1; %Will perturb D-matrix
  [VE,eigE]=eig(SERE); invVE=VE^(-1); eigE=diag(eigE);
else
  Eflag=0;  
end

SERCnew=SERC; SERDnew=SERD; SEREnew=SERE;


N=length(SERA);

bigB=[];
bigc=[];
Ns=length(s);
Ns2=length(s2);
Nc=length(SERD);
Nc2=Nc*Nc;
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



    
if isempty(RPopts.H) 

if RPopts.outputlevel==1, disp('    Building system equation (once)...'),end 
%========================================
% LOOP FOR LEAST SQUARES PROBLEM:
%========================================
if (Dflag+Eflag)==2
  bigA=zeros(Ns*Nc2,Nc*(N+2));    
elseif (Dflag+Eflag)==1    
  bigA=zeros(Ns*Nc2,Nc*(N+1));
else    
  bigA=zeros(Ns*Nc2,Nc*(N));
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


for k=1:Ns
      sk=s(k);
      % Calculating matrix Mmat (coefficients for LS-problem)
      tell=0;
      offs=0;
      Yfit=fitcalcPRE(sk,SERA,SERC,SERD,SERE);
      %y=SERD(row,col)+sk*SERE(row,col);
      %y=y+sum(squeeze(SERC(row,col,1:N))./(s(k)-SERA(1:N)));if
      
      
      if ~isempty(weightparam)
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
      else
        weight=bigweight(:,:,k);  
      end    
      
  for m=1:N
    V=squeeze(bigV(:,(m-1)*Nc+1:m*Nc));
    invV=V^(-1);%squeeze(biginvV(:,(m-1)*Nc+1:m*Nc)); 
    if cindex(m)==0 %real pole
      dum=1/(sk-SERA(m));   
    elseif cindex(m)==1  %complex pole, 1st part
      dum=(1/(sk-SERA(m)) + 1/(sk-SERA(m)') );   
    else
      dum=(i/(sk-SERA(m)') - i/(sk-SERA(m)) ); 
    end 
    
    for egenverdi=1:Nc
      tell=0;
      gamm=V(:,egenverdi)*invV(egenverdi,:);     
      for row=1:Nc
        %ind2=0;
        for col=1:Nc
          faktor=weight(row,col);    
          tell=tell+1;
          if cindex(m)==0 %real pole
            Mmat(tell,offs+egenverdi)=gamm(row,col)*faktor*dum;   
          elseif cindex(m)==1  %complex pole, 1st part
            Mmat(tell,offs+egenverdi)=gamm(row,col)*faktor*dum;   
          else
            Mmat(tell,offs+egenverdi)=gamm(row,col)*faktor*dum; 
          end    
        end  
      end  
    end %for egenverdi=1:Nc
    offs=offs + (Nc);
  end %for m=1:N
 
  if Dflag==1
    for egenverdi=1:Nc %Egenverdi til residue-matrise
      gamm=VD(:,egenverdi)*invVD(egenverdi,:); %Ytre-produkt 
      tell=0;
      for row=1:Nc
        for col=1:Nc
          tell=tell+1;  
          faktor=weight(row,col);  
          Mmat(tell,offs+egenverdi)=gamm(row,col)*faktor; 
        end
      end  
    end  
  end  
  
  if Eflag==1
    for egenverdi=1:Nc %Egenverdi til residue-matrise
      gamm=VE(:,egenverdi)*invVE(egenverdi,:); %Ytre-produkt 
      tell=0;
      for row=1:Nc
        for col=1:Nc
          tell=tell+1;  
          faktor=weight(row,col);  
          Mmat(tell,offs+Nc*Dflag+egenverdi)=gamm(row,col)*sk*faktor; 
        end
      end  
    end  
  end  
  bigA((k-1)*Nc2+1:k*Nc2,:)=Mmat;  
 
  
  
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


%bigA2=zeros(Ns4*Nc2,Nc*(N));
bigA2=zeros(Ns4*Nc2,Nc*(N+Dflag+Eflag));
 
for k=1:1:Ns4
      sk=s4(k);
      % Calculating matrix Mmat (coefficients for LS-problem)
      tell=0;
      offs=0;
      Yfit=fitcalcPRE(sk,SERA,SERC,SERD,SERE);
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
      dum=1/(sk-SERA(m));   
    elseif cindex(m)==1  %complex pole, 1st part
      dum=(1/(sk-SERA(m)) + 1/(sk-SERA(m)') );   
    else
      dum=(i/(sk-SERA(m)') - i/(sk-SERA(m)) ); 
    end 
    
    for egenverdi=1:Nc
      tell=0;
      gamm=V(:,egenverdi)*invV(egenverdi,:);     
      for row=1:Nc
        %ind2=0;
        for col=1:Nc
          faktor=weight(row,col);    
          tell=tell+1;
          if cindex(m)==0 %real pole
            Mmat(tell,offs+egenverdi)=gamm(row,col)*faktor*dum;   
          elseif cindex(m)==1  %complex pole, 1st part
            Mmat(tell,offs+egenverdi)=gamm(row,col)*faktor*dum;   
          else
            Mmat(tell,offs+egenverdi)=gamm(row,col)*faktor*dum; 
          end    
        end  
      end  
    end %for egenverdi=1:Nc
    offs=offs + (Nc);
  end %for m=1:N
 
  if Dflag==1
    for egenverdi=1:Nc %Egenverdi til residue-matrise
      gamm=VD(:,egenverdi)*invVD(egenverdi,:); %Ytre-produkt 
      tell=0;
      for row=1:Nc
        for col=1:Nc
          tell=tell+1;  
          faktor=weight(row,col);  
          Mmat(tell,offs+egenverdi)=gamm(row,col)*faktor; 
        end
      end  
    end  
  end  
  
  if Eflag==1
    for egenverdi=1:Nc %Egenverdi til residue-matrise
      gamm=VE(:,egenverdi)*invVE(egenverdi,:); %Ytre-produkt 
      tell=0;
      for row=1:Nc
        for col=1:Nc
          tell=tell+1;  
          faktor=weight(row,col);  
          Mmat(tell,offs+Nc*Dflag+egenverdi)=gamm(row,col)*sk*faktor; 
        end
      end  
    end  
  end  
  bigA2((k-1)*Nc2+1:k*Nc2,:)=Mmat; 
  
end %for k=1:Ns
bigA=[bigA;bigA2];

end %if auxflag==1


bigA=[real(bigA); imag(bigA)];
Acol=length(bigA(1,:));
for col=1:Acol
  Escale(col)=norm(bigA(:,col),2);
  bigA(:,col)=bigA(:,col)./Escale(col);
end
H=bigA.'*bigA;
RPopts.H=H; clear H
RPopts.Escale=Escale; clear Escale
RPopts.bigV=bigV;
RPopts.biginvV=biginvV;
if RPopts.outputlevel==1, disp('    Done'), end 

else %if isempty(RPopts)
  bigV=RPopts.bigV;  
  biginvV=RPopts.biginvV;
  if Dflag~=RPopts.oldDflag | Eflag~=RPopts.oldEflag
    RPopts.H=RPopts.H(1:Nc*(N+Dflag+Eflag),1:Nc*(N+Dflag+Eflag));
    RPopts.Escale=RPopts.Escale(1:Nc*(N+Dflag+Eflag));
  end  
end %if isempty(RPopts) 



Mmat2=zeros(Nc2,Nc*(N+Dflag+Eflag));
viol_G=[];
viol_D=[];
viol_E=[];
%========================================
% LOOP FOR CONSTRAINT PROBLEM, TYPE #1 (violating eigenvalues in s2):
%========================================
for k=1:Ns2
  sk=s2(k);
  for row=1:Nc
    for col=1:Nc
      Y(row,col)=SERD(row,col)+sk*SERE(row,col);
      Y(row,col)=Y(row,col)+sum(squeeze(SERC(row,col,1:N))./(sk-SERA(1:N)));
    end
  end 

  %Calculating eigenvalues and eigenvectors: 
  [V,Z] = eig(real(Y)); Z=diag(Z);
  EE(:,k)=real(Z);
  if min(real(Z))<0  %any violations

    % Calculating matrix M2mat; matrix of partial derivatives :
    tell=0;
    offs=0;


  for m=1:N
    VV=bigV(:,(m-1)*Nc+1:m*Nc);
    invVV=biginvV(:,(m-1)*Nc+1:m*Nc);       
    for egenverdi=1:Nc
      tell=0;
      gamm=VV(:,egenverdi)*invVV(egenverdi,:);     
      for row=1:Nc
        %ind2=0;
        for col=1:Nc
          tell=tell+1;
          if cindex(m)==0 %real pole
            Mmat2(tell,offs+egenverdi)=gamm(row,col)/(sk-SERA(m));   
          elseif cindex(m)==1  %complex pole, 1st part
            Mmat2(tell,offs+egenverdi)=gamm(row,col)*(1/(sk-SERA(m)) + 1/(sk-SERA(m)') );   
          else
            Mmat2(tell,offs+egenverdi)=gamm(row,col)*(i/(sk-SERA(m)') - i/(sk-SERA(m)) ); 
          end    
        end  
      end  
    end %for egenverdi=1:Nc
    offs=offs + (Nc);
  end %for m=1:N
 
  if Dflag==1
    for egenverdi=1:Nc %Egenverdi til residue-matrise
      tell=0;  
      gamm=VD(:,egenverdi)*invVD(egenverdi,:); %Ytre-produkt
      for row=1:Nc
        for col=1:Nc  
          tell=tell+1;
          Mmat2(tell,offs+egenverdi)=gamm(row,col);
        end  
      end
    end  
  end 
    
  
    for n=1:Nc
      tell=0;
      V1=V(:,n); 
      for row=1:Nc
        for col=1:Nc
          if row==col
            qij=V1(row)^2;
          else
            qij=V1(row)*V1(col); 
          end
          tell=tell+1;
          Q(n,tell)=qij;
        end
      end
    end


    B=Q*Mmat2; 
    delz=real(Z);
    for n=1:Nc   %ustabilitet?
      if delz(n)<0
        bigB=[bigB;B(n,:)];
        bigc=[bigc;-TOL+delz(n)];
        viol_G=[viol_G; delz(n)];
      end
    end
          
  end %if max(real(Z))>0  %any violations
  
end %k



%========================================
% LOOP FOR CONSTRAINT PROBLEM, TYPE #2: (all eigenvalues in s3):
%========================================
Ns3=length(s3);
for k=1:Ns3
  sk=s3(k);
  for row=1:Nc
    for col=1:Nc
      Y(row,col)=SERD(row,col)+sk*SERE(row,col);
      Y(row,col)=Y(row,col)+sum(squeeze(SERC(row,col,1:N))./(sk-SERA(1:N)));
    end
  end 

  %Calculating eigenvalues and eigenvectors: 
  [V,Z] = eig(real(Y)); Z=diag(Z);
  EE(:,k)=real(Z);
  %if min(real(Z))<0  %any violations

    % Calculating matrix M2mat; matrix of partial derivatives :
    tell=0;
    offs=0;


  for m=1:N
    VV=bigV(:,(m-1)*Nc+1:m*Nc);
    invVV=biginvV(:,(m-1)*Nc+1:m*Nc);       
    for egenverdi=1:Nc
      tell=0;
      gamm=VV(:,egenverdi)*invVV(egenverdi,:);     
      for row=1:Nc
        %ind2=0;
        for col=1:Nc
          tell=tell+1;
          if cindex(m)==0 %real pole
            Mmat2(tell,offs+egenverdi)=gamm(row,col)/(sk-SERA(m));   
          elseif cindex(m)==1  %complex pole, 1st part
            Mmat2(tell,offs+egenverdi)=gamm(row,col)*(1/(sk-SERA(m)) + 1/(sk-SERA(m)') );   
          else
            Mmat2(tell,offs+egenverdi)=gamm(row,col)*(i/(sk-SERA(m)') - i/(sk-SERA(m)) ); 
          end    
        end  
      end  
    end %for egenverdi=1:Nc
    offs=offs + (Nc);
  end %for m=1:N
  if Dflag==1
    for egenverdi=1:Nc %Egenverdi til residue-matrise
      tell=0;  
      gamm=VD(:,egenverdi)*invVD(egenverdi,:); %Ytre-produkt
      for row=1:Nc
        for col=1:Nc  
          tell=tell+1;
          Mmat2(tell,offs+egenverdi)=gamm(row,col);
        end  
      end
    end  
  end 
  
    for n=1:Nc
      tell=0;
      V1=V(:,n); 
      for row=1:Nc
        for col=1:Nc
          if row==col
            qij=V1(row)^2;
          else
            qij=V1(row)*V1(col); 
          end
          tell=tell+1;
          Q(n,tell)=qij;
        end
      end
    end

    B=Q*Mmat2; 
    delz=real(Z);
    for n=1:Nc   
      bigB=[bigB;B(n,:)];
      bigc=[bigc;-TOL+delz(n)];
      viol_G=[viol_G; delz(n)];     
    end
          
  %end %if max(real(Z))>0  %any violations
  
end %for k=1:Ns3

%Legger til constraint for evt. e-verdier i D <0
if Dflag==1
  for n=1:Nc
    %if eigD(n)<0
      dum=zeros(1,(Nc*(N+Dflag+Eflag)));
      dum(Nc*N+n)=1;
      bigB=[bigB;dum];
      bigc=[bigc; eigD(n)-TOL];
      viol_G=[viol_G; eigD(n)];      
      viol_D=[viol_D; eigD(n)];          
    %end
  end
end

%Legger til constraint for evt. e-verdier i E <0
if Eflag==1
  for n=1:Nc
    %if eigD(n)<0
      dum=zeros(1,(Nc*(N+Dflag+Eflag)));
      dum( Nc*(N+Dflag)+n )=1;
      bigB=[bigB;dum];
      bigc=[bigc; eigE(n)-TOLE];
      viol_E=[viol_E; eigE(n)];          
    %end
  end
end



if length(bigB)==0 
  return   %No passivity violations
end


c=bigc;



bigB=[real(bigB)];
for col=1:length(RPopts.H)
  if length(bigB)>0
    bigB(:,col)=bigB(:,col)./RPopts.Escale(col);
  end
end  



ff=zeros(length(RPopts.H),1);

clear bigA;
%disp('Starting quadprog...')

warning off


if RPopts.solver=='QUADPROG'
  opts.Display='off';
  [dx,lambda]=quadprog(RPopts.H,ff,-bigB,bigc,[],[],[],[],[],opts);
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



dx=dx./RPopts.Escale.';


%Oppdaterer egenverdiene og SERC:
for m=1:N
  if cindex(m)==0 %real pole
    D1=diag(dx((m-1)*Nc+1:m*Nc));  
    SERCnew(:,:,m)=SERCnew(:,:,m) +bigV(:,(m-1)*Nc+1:m*Nc)*D1*biginvV(:,(m-1)*Nc+1:m*Nc);  
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
if Eflag==1
  EE=diag(dx((N+Dflag)*Nc+1:(N+Dflag+Eflag)*Nc));
  SEREnew=SEREnew+VE*EE*invVE;
end

%==========================
SERDnew=(SERDnew+SERDnew.')/2;
SEREnew=(SEREnew+SEREnew.')/2;
for m=1:N
  SERCnew(:,:,m)=(SERCnew(:,:,m)+SERCnew(:,:,m).')/2;
end  
%==========================
SER.R=SERCnew;
SER.D=SERDnew;
SER.E=SEREnew;
[SER]=pr2ss(SER);  

RPopts.oldDflag=Dflag;
RPopts.oldEflag=Eflag;

end %function FRPY 


