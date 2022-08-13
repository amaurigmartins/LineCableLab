function [SER1,bigYfit,opts2]=RPdriver(SER,s,opts)
% SYNTAX
%    [SERnew,opts2]=RPdriver(SER,s)      %use default settings 
%    [SERnew,opts2]=RPdriver(SER,s,opts) %include structure "opts" to override
%    default settings
%
%      ======================================================
%      =   Routine: RPdriver.m                              = 
%      =   Version 1.0                                      =
%      =   Last revised: 06.10.2008                         = 
%      =   Bjorn Gustavsen                                  =
%      =   SINTEF Energy Research, N-7465 Trondheim, NORWAY =
%      =   bjorn.gustavsen@sintef.no                        =
%      ======================================================
%
%  function [SERnew]=RPdriver(SER,s,opts);
%  
%  PURPOSE : Perturb eigenvalues of R, D, and E of an Nth order pole-residue model with Nc ports
%                  
%                   N   Rm
%           Y(s)=SUM(---------  ) +D +s*E 
%                  m=1 (s-am)        
%
%            - to enforce passivity: eig(real(Y))>0 for all frequencies
%            - to enforce a positive definite D (asymptotic passivity)
%            - to enforce a positive definite E 
%
%  This ensures that the model will not cause instabilities when used in a time domain 
%  simulation. 
%
% =================================================
%     INPUT :
% =================================================
%
%  SER.poles(1,N)     : Poles
%  SER.R    (Nc,Nc,N) : Residues 
%  SER.D    (Nc,Nc)   : Constant term
%  SER.E    (Nc,Nc)   : Capacitive term (y-parameter model) 
%  s        (1,Ns)    : Frequency samples (jw [rad/sec])
%
%
% The following additional input parameters can be specified to override default settings: 
%
% opts.parametertype='Y' --> Assume Y-parameter model 
% opts.parametertype='S' --> Assume S-parameter model 
%   Default value: 'Y'
%
% opts.Niter_out : Maximum n.o. iterions in outer-loop of "Robust Iterations" (Default: 10)
% opts.Niter_in  : Maximum n.o. iterions in inner-loop of "Robust Iterations" (Default: 0)
% opts.TOLGD     : Target value for negative eigenvalues of G and D (Default: 1E-6)     
% opts.TOLE      : Target value for negative eigenvalues of E       (Default: 1E-12)  
%
% opts.cmplx_ss=1 --> Will return a complex state-space model with diagonal A  (Def=1)
% opts.cmplx_ss=0 --> Will return a real-only state-space model with block-diagonal A
%
% opts.weightparam (Def=1)
%            =1 --> weight=1 for all elements in Least Sq. problem, at all freq. 
%            =2 --> weight(s)=1/abs(Yij(s))      ; indvidual element weight
%            =3 --> weight(s)=1/sqrt(abs(Yij(s))); indvidual element weight
%            =4 --> weight(s)=1/norm(Y(s))       ; common weight for all matrix elements
%            =5 --> weight(s)=1/sqrt(norm(Y(s)) ; common weight for all matrix elements
%
% opts.weightfactor: Weighting for auxiliary samples at out-of-band frequencies (pole-locations).
%                    Improves numerical robustness. Default value: 1e-3
%
% opts.colinterch  : =1 --> Will retain the correct sequence of eigenvalues (Y-paramerters) or 
%                           singular values (S-parameters) as function of frequency during assessment 
%                           of passivity violation extrema. (Def=1)
%
% opts.outputlevel : =0 --> Minimum level of output information to screen
%                    =1 --> Maximum level of output information to screen (Def=1)
%
% opts.plot        : structure requesting to plot eigenvalues of G(s)=real(Y(s)), during iterations. 
% opts.plot.s_pass : row-vector containing frequency samples (s=jw) for plotting eigenvalues 
% opts.plot.xlim   : xlimits for plots, in Hz. (e.g. [1 1000])
% opts.plot.ylim   : ylimits for plots (e.g. [-0.1 0.1])
%
% =================================================
%     OUTPUT :
% =================================================
% - SER1 : New structure with perturbed SER.R, SER.D, and SER.E, as well as SER.A and SER.C
%          of the associated state-space model.
% - opts2: Structure containing all settings used in the call (including default ones)
% =================================================
%     APPROACH :
% =================================================
% The perturbation is done using the "Fast Residue Perturbation" technique 
% described in [1],[3], combined with passivity checking via a half-size test
% matrix [2],[4].
%
%***************************************************************************** 
% NOTE: The use of this program is limited to NON-COMMERCIAL usage only.
% If the program code (or a modified version) is used in a scientific work, 
% then two references should be made, depending on the application:  
%
% a) When applied to Y-parameter data sets:
%                      
% [1] B. Gustavsen, "Fast passivity enforcement for pole-residue models by perturbation
%     of residue matrix eigenvalues", 
%     IEEE Trans. Power Delivery, vol. 23, no. 4, pp. 2278-2285, October 2008. 
% 
% [2] A. Semlyen and B. Gustavsen, "A half-size singularity test matrix for fast       
%     and reliable passivity assessment of rational models",        
%     IEEE Trans. Power Delivery, accepted for publication.      
%
% b) When applied to S-parameter data sets:
%
% [3] B. Gustavsen, "Fast passivity enforcement for S-parameter models by perturbation
%     of residue matrix eigenvalues", 
%     IEEE Trans. Advanced Packaging, accepted for publication. 
%
% [4] B. Gustavsen and A. Semlyen, "Fast passivity assessment for S-parameter rational
%     models via a half-size test matrix", 
%     IEEE Trans. Microwave Theory and Techniques, accepted for publication. 
%
%**************************************************************************
%
% This file is part of the Matrix Fitting Toolbox, v1.
% Filename: RPdriver.m
% Package: Matrix_Fitting_Toolbox_1.zip.


%--------------------------------------------
MPopts.auxflag     =1;
MPopts.solver      ='QUADPROG';
%--------------------------------------------- 

[SER]=pr2ss(SER);  %Convert model from pole-residue to state-space

    
disp('-----------------S T A R T--------------------------')

if nargin<3,opts=[];end
    

if ~isfield(opts,'parametertype'),   opts.parametertype='Y'; end  
if ~isfield(opts,'Niter_out'),       opts.Niter_out=10; end  
if ~isfield(opts,'Niter_in'),        opts.Niter_in=0; end  
if ~isfield(opts,'TOLGD'),           opts.TOLGD=1e-6; end  
if ~isfield(opts,'TOLE'),            opts.TOLE=1e-12; end  
if ~isfield(opts,'cmplx_ss'),        opts.cmplx_ss=1; end  
if ~isfield(opts,'weightfactor'),    opts.weightfactor=0.001; end 
if ~isfield(opts,'weightparam'),     opts.weightparam=1; end  
if ~isfield(opts,'method'),          opts.method='FRP'; end  
if ~isfield(opts,'colinterch'),      opts.colinterch=1; end
if ~isfield(opts,'outputlevel'),     opts.outputlevel=1; end 
if ~isfield(opts,'weight'),     opts.weight=[]; end 


if opts.method=='FMP' & opts.parametertype=='S'
  %disp('ERROR in RPdriver.m: FMP cannot be used together with S-parameters. Must stop.')
  %return
end  
colinterch=opts.colinterch;

MPopts.TOLGD=opts.TOLGD;
MPopts.TOLE=opts.TOLE;
MPopts.weightfactor=opts.weightfactor;
MPopts.weightparam=opts.weightparam;
MPopts.weight=opts.weight;
MPopts.outputlevel=opts.outputlevel;

opts2=struct('method',opts.method,'parametertype',opts.parametertype,'Niter_out',opts.Niter_out,'Niter_in',opts.Niter_in,...
             'TOLGD',opts.TOLGD,'TOLE',opts.TOLE,'cmplx_ss',opts.cmplx_ss,'weightparam',opts.weightparam,...
             'weightfactor',opts.weightfactor,'colinterch',opts.colinterch,'outputlevel',opts.outputlevel);


if opts.parametertype=='Y'
  disp('*** Y-PARAMETERS ***')    
elseif opts.parametertype=='S'
  disp('*** S-PARAMETERS ***')    
end    
         
plotte=0;
if isfield(opts, 'plot')
  if ~isempty(opts.plot)  
    plotte=1;
    s_pass=opts.plot.s_pass; %For plotting purposes
    if isfield(opts.plot,'xlim'),xlimflag=1;else,xlimflag=0;end
    if isfield(opts.plot,'ylim'),ylimflag=1;else,ylimflag=0;end  
  end
end
break_outer=0;
olds3=[];

SER0=SER;
Nc=length(SER.D);

Niter_out=opts.Niter_out;
Niter_in=opts.Niter_in;

%========================================================
%Plotting eigenvalues of original model (SERC0, SERD0): 
%========================================================
if plotte==1
  oldT0=[];
  oldU=[];
  I=sparse(ones(length(SER.A(:,1)),1));
  for k=1:length(s_pass)
    Y   =SER.C*diag( ( s_pass(k)*I -diag(SER.A) ).^(-1) )*SER.B +SER.D +s_pass(k)*SER.E;
    if opts.parametertype=='Y'
      G=real(Y); 
      [T0,D]=eig(G);
      T0=rot(T0); %minimizing phase angle of eigenvectors in least squares sense  
      [T0,D]=intercheig(T0,oldT0,D,Nc,k); oldT0=T0;
      EE0(:,k)=diag(D);
    elseif opts.parametertype=='S'
      if colinterch==0  
        EE0(:,k)=svd(Y,0); 
      else 
        [U,S,V]=svd(Y,0);   
        [U,S,V]=interchsvd(U,oldU,S,V,Nc,k); oldU=U;   
        EE0(:,k)=diag(S);   
      end  
    end    
  end
  figure(7),h0=plot(s_pass/(2*pi*i),EE0.','b'); 
  %figure(7),h0=plot(s_pass/(2*pi*i),EE0.','k');   
  hold on, grid on,
  if xlimflag==1, xlim([opts.plot.xlim]); else, xlim([s_pass(1)/(2*pi*i) s_pass(end)/(2*pi*i)]); end
  if ylimflag==1, ylim([opts.plot.ylim]);end
  xlabel('Frequency [Hz]'); 
  if opts.parametertype=='Y'
    ylabel('Eigenvalues of G_{rat}')
  else  
    ylabel('Singular values of S_{rat}')  
  end  
  %!!title('Eigenvalues of G(s)');

  figure(8),h2=plot(s_pass/(2*pi*i),EE0.','b');  
  hold on, grid on, 
  if xlimflag==1, xlim([opts.plot.xlim]); else, xlim([s_pass(1)/(2*pi*i) s_pass(end)/(2*pi*i)]); end
  if ylimflag==1, ylim([opts.plot.ylim]);end
  title('Monitoring enforcement process (eig(G(s))'); xlabel('Frequency [Hz]'); 
  grid on; pause(0.1); 
end %if plotte==1


outputlevel=opts.outputlevel;
%if outputlevel==1, %Stopwatch only with high-level output
%  tic;
%end  
t=[0 0 0 0];
%t(1): assessing crossover frequencies 
%t(2): assessing violation extrema  
%t(3): passivity enforcement
%t(4): total

% tdumtotal=toc;



%=======================================================
% Passivity Enforcement :
%=======================================================
QP.first=1;
QPopts=[];

SER1=SER0;

for iter_out=1:Niter_out
 
  if break_outer==1
    SER0=SER1;  
    break
  end
  s3=[];  %Initial value  
      
      
  for iter_in=1:1:Niter_in+1
  
  s2=[];  
  SERflag=1;
   if outputlevel==1 
     disp(' ')
     disp([     '  [ ' num2str(iter_out) '  ' num2str(iter_in-1) ' ]']) 
     disp( '  Passivity Assessment:')   
   end
%    tdum=toc;
   
   if iter_in==1, %######
   if opts.parametertype=='Y'
     [wintervals]=pass_check_Y(SERflag,SER.poles,[],SER1.R,SER1.D,colinterch);
   else  
     TOL=1e-3; spy=1;  
     [wintervals]=pass_check_S(SERflag,SER.poles,[],SER1.R,SER1.D,TOL,spy,colinterch);  
   end  
%    t(1)=t(1) +toc-tdum;
  
    if length(wintervals)>0  
      if outputlevel==1  
        disp(['    N.o. violating intervals: ' num2str(length(wintervals(1,:)))]) 
      end
    end
   

  if opts.parametertype=='Y'
    if length(wintervals)==0 & all(eig(SER1.D)>=0) & all(eig(SER1.E)>=0); SER0=SER1; break_outer=1; break, end
  elseif opts.parametertype=='S'     
    if length(wintervals)==0 & all(svd(SER1.D)<=1); SER0=SER1; break_outer=1; break, end
  end   

  if opts.parametertype=='Y'
    %Identifying minima within each interval:
%     tdum=toc;  
    [s_viol,g_pass,ss]=violextremaY(SERflag,wintervals.',SER.poles,[],SER1.R,SER1.D,colinterch);
%     t(2)=t(2) +toc-tdum; 
    s2=[s_viol.'];s2=sort(s2);
    if length(s2)==0 & all(eig(SER1.D)>0)  
      break 
    end
  elseif opts.parametertype=='S'  
     %Identifying maxima within each interval:
%     tdum=toc;  
    [s_viol,g_pass,ss]=violextremaS(SERflag,wintervals.',SER.poles,[],SER1.R,SER1.D,colinterch);
%     t(2)=t(2) +toc-tdum;
    s2=[s_viol.'];s2=sort(s2);
    SER1=SER0; 
    if length(s2)==0 & all(svd(SER1.D)<1) 
      break 
    end
  end    
%s2=2*pi*i*1e2;
 
   end %if iter_in==1, %######


  
    if iter_in==1 
        
      if outputlevel==1  
      if opts.parametertype=='Y' 
        if min(g_pass)<0   
          disp(['    Max. violation, eig(G) : ' num2str(g_pass) ' @ ' num2str((ss)/(2*pi*i)) ' Hz'])
        else
          disp(['    Max. violation, eig(G) :  None']) 
        end
        if min(eig(SER0.D))<0   
          disp(['    Max. violation, eig(D) : ' num2str(min(eig(SER1.D)))])
        else
          disp(['    Max. violation, eig(D) :  None']) 
        end
         if min(eig(SER0.E))<0   
          disp(['    Max. violation, eig(E) : ' num2str(min(eig(SER1.E)))])
        else
          disp(['    Max. violation, eig(E) :  None']) 
        end        
      
      elseif opts.parametertype=='S'
        if max(g_pass)>1   
          disp(['    Max. violation, sing(S) : ' num2str(-1+g_pass) ' @ ' num2str((ss)/(2*pi*i)) ' Hz'])
        else
          disp(['    Max. violation, sing(S) :  None']) 
        end
        if max(svd(SER0.D))>1   
          disp(['    Max. violation, sing(D) : ' num2str(-1+max(svd(SER1.D)))])
        else
          disp(['    Max. violation, sing(D) :  None']) 
        end      
      end
      end %if outputlevel==1
      if outputlevel~=1
        if opts.parametertype=='Y' 
          min1=min(g_pass);min2=min(eig(SER1.D));  
          disp(['    Max. violation  : ' num2str(min([min1 min2]))])
          if min(eig(SER0.E))<0
          disp(['    Max. violation, E: ' num2str(min(eig(SER1.E)))])
          end
        elseif opts.parametertype=='S' 
          max1=max(g_pass);max2=max(eig(SER1.D));  
          disp(['    Max. violation  : ' num2str(-1+max([max1 max2]))])            
        end    
      end    
    end %if iter_in==1     
     
    if length(s3)>0
      %!!SER1=SER0; 
    end
    %disp('Starting Passivity enforcement...')

   %[SER1,MPopts]=FMP(SER1,s,s2,s3,MPopts); 
   if outputlevel==1, disp( '  Passivity Enforcement...'), end
%    tdum=toc;
   if opts.method=='FMP'     
     [SER1,MPopts]=FMP(SER0,s,s2,s3,MPopts);  
   elseif opts.method=='FRP'
     if opts.parametertype=='Y'
       [SER1,MPopts]=FRPY(SER0,s,s2,s3,MPopts);    
     else    
       [SER1,MPopts]=FRPS(SER0,s,s2,s3,MPopts);
     end  
   else
     disp('****** ERROR #1 in FMPdriver.m')  
   end  
%    t(3)=t(3) +toc-tdum;

%     t(4)=t(4) +toc-tdumtotal;
    if plotte==1
      if opts.parametertype=='Y' 
        oldT0=[];
        tell=0;
        I=sparse(ones(length(SER.A(:,1)),1));
        for k=1:length(s_pass)
          Y   =SER1.C*diag( ( s_pass(k)*I -diag(SER1.A) ).^(-1) )*SER1.B +SER1.D +s_pass(k)*SER1.E; 
          G=real(Y); 
          [T0,D]=eig(G);
          T0=rot(T0); %minimizing phase angle of eigenvectors in least squares sense  
          [T0,D]=intercheig(T0,oldT0,D,Nc,k); oldT0=T0;
          EE1(:,k)=diag(D);
        end
      elseif opts.parametertype=='S'
        oldU=[];
        tell=0;
        I=sparse(ones(length(SER.A(:,1)),1));
        for k=1:length(s_pass)
          Y   =SER1.C*diag( ( s_pass(k)*I -diag(SER1.A) ).^(-1) )*SER1.B +SER1.D +s_pass(k)*SER1.E;
          if colinterch==0  
            EE1(:,k)=svd(Y,0); 
          else 
            [U,S,V]=svd(Y,0);   
            [U,S,V]=interchsvd(U,oldU,S,V,Nc,k); oldU=U;   
            EE1(:,k)=diag(S);   
          end 
        end    
      end    
      figure(8),h2=plot(s_pass/(2*pi*i),EE0.','b-'); hold on
      figure(8),h3=plot(s_pass/(2*pi*i),EE1.','r--');
      if xlimflag==1, xlim([opts.plot.xlim]); else, xlim([s_pass(1)/(2*pi*i) s_pass(end)/(2*pi*i)]); end
      if ylimflag==1, ylim([opts.plot.ylim]);end
      title('Monitoring enforcement process'); xlabel('Frequency [Hz]');
      hold off; grid on, 
      legend([h2(1) h3(1)],'Previous','Perturbed'); drawnow, pause(0.01); hold off
    end
    
%     tdumtotal=toc;
    if iter_in~=Niter_in+1 %Not last run in inner-loop
      if opts.parametertype=='Y'  
%         tdum=toc;
        [wintervals]=pass_check_Y(SERflag,SER1.poles,[],SER1.R,SER1.D);
%         t(1)=t(1) +toc-tdum;   
%         tdum=toc;      
        [s_viol]=violextremaY(SERflag,wintervals.',SER1.poles,[],SER1.R,SER1.D,colinterch);  
%         t(2)=t(2) +toc-tdum;        
      elseif opts.parametertype=='S'  
%          tdum=toc;  
         [wintervals]=pass_check_S(SERflag,SER1.poles,[],SER1.R,SER1.D,TOL,spy,colinterch);
%           t(1)=t(1) +toc-tdum;        
%          tdum=toc;  
         [s_viol]=violextremaS(SERflag,wintervals.',SER1.poles,[],SER1.R,SER1.D,colinterch);  
%          t(2)=t(2) +toc-tdum;
      end  
      olds3=s3;
      s3=[s3; s2; s_viol.'];
    end
      
   
    %!!if length(s3)==length(olds3) | iter_in==Niter_in+1 %No new violations or last iteration
    if iter_in==Niter_in+1 %last inner-iteration
      s3=[]; s2=[]; 
      if plotte==1
        EE0=EE1; %Update model
      end  
      SER0=SER1;
      %break %leave inner loop
    end 

%     t(4)=t(4) +toc-tdumtotal;    

  end %for iter=1:1


  
end %for iter_out=1:2


%========================================================
%Plotting eigenvalues of modified model (SERC1, SERD1): 
%=========================================================
if plotte==1
  if opts.parametertype=='Y'     
    oldT0=[];
    tell=0;
    for k=1:length(s_pass)
      Y   =SER1.C*diag( ( s_pass(k)*I -diag(SER1.A) ).^(-1) )*SER1.B +SER1.D +s_pass(k)*SER1.E;  
      G=real(Y); 
      [T0,D]=eig(G);
      T0=rot(T0); %minimizing phase angle of eigenvectors in least squares sense  
      [T0,D]=intercheig(T0,oldT0,D,Nc,k); oldT0=T0;
      EE1(:,k)=diag(D);
    end
  elseif opts.parametertype=='S'
    oldU=[];
    tell=0;
    for k=1:length(s_pass)
      Y   =SER1.C*diag( ( s_pass(k)*I -diag(SER1.A) ).^(-1) )*SER1.B +SER1.D +s_pass(k)*SER1.E;  
      if colinterch==0  
        EE1(:,k)=svd(Y,0); 
      else 
        [U,S,V]=svd(Y,0);   
        [U,S,V]=interchsvd(U,oldU,S,V,Nc,k); oldU=U;   
        EE1(:,k)=diag(S);   
      end 
    end    
  end    
      
  figure(7),h1=plot(s_pass/(2*pi*i),EE1.','r--');  
  %a=2; figure(7),h1=plot(s_pass(1:a:end)/(2*pi*i),EE1(:,1:a:end).','k.');
  if xlimflag==1, xlim([opts.plot.xlim]); else, xlim([s_pass(1)/(2*pi*i) s_pass(end)/(2*pi*i)]); end
  if ylimflag==1, ylim([opts.plot.ylim]);end
  legend([h0(1) h1(1)],'Original','Perturbed',1); hold off
  set(gcf,'paperunits','centimeters');
  %newpos =[0.25 2.5 28 13]; %NB: I Word: resize height til 80%.
  newpos =[0.25 2.5 9.5 9];  %NB: I Word: resize height til 80%. 
  set(gcf,'paperpos',newpos);

end %if plotte==1  


if isempty(wintervals)
  if outputlevel==1  
  disp(' ')
  end
  disp('-->Passivity was successfully enforced.')  
  if outputlevel==1
  if opts.parametertype=='Y'
    disp(['   Max. violation, eig(G) :  None'])
    disp(['   Max. violation, eig(D) :  None'])    
    disp(['   Max. violation, eig(E) :  None'])  
  elseif opts.parametertype=='S'
    disp(['   Max. violation, sing(S) :  None'])
    disp(['   Max. violation, sing(D) :  None'])    
  end
  end
else
    %if iter_in==1 
      disp(['   ***Max. violation, eig(G) : ' num2str(min(g_pass))])
      disp(['   ***Max. violation, eig(D) : ' num2str(min(eig(SER0.D)))])
      disp(['   ***Max. violation, eig(E) : ' num2str(min(eig(SER0.E)))])      
    %end    
  disp('-->Iterations terminated before completing passivity enforcement.')
  disp('   Increase parameter opts.Niter_out.')  
end    

%Producing plot;
Ns=length(s);
bigYfit=zeros(Nc,Nc,Ns);
I=sparse(ones(length(SER.A(:,1)),1));
for k=1:Ns
  Y =SER1.C*diag( ( s(k)*I -diag(SER1.A) ).^(-1) )*SER1.B +SER1.D +s(k)*SER1.E;
  bigYfit(:,:,k)=Y;
end

%Converting to real-only state-space, if requested
if opts.cmplx_ss==0
  N=length(SER1.A);
  %SER1.A=diag(sparse(SER1.A));
  cindex=zeros(1,N);
  for m=1:N 
    if imag(SER1.A(m,m))~=0  
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
  n=0;
  for m=1:N
    n=n+1;
    if cindex(m)==1
      a=SER1.A(n,n); a1=real(a); a2=imag(a);
      c=SER1.C(:,n); c1=real(c); c2=imag(c);
      b=SER1.B(n,:); b1=2*real(b); b2=-2*imag(b);
      Ablock=[a1 a2;-a2 a1];
   
      SER1.A(n:n+1,n:n+1)=Ablock;
      SER1.C(:,n)=c1;
      SER1.C(:,n+1)=c2;
      SER1.B(n,:)=b1;
      SER1.B(n+1,:)=b2;
    end
  end   
end %if opts.opts.cmplx_ss==0

dum=t(1)+t(2);
if outputlevel==1
disp(['Time summary: '])
disp(['   Passivity assessment : ' num2str(dum) ' sec'])
disp(['   Passivity enforcement: ' num2str(t(3)) ' sec'])
disp(['   Total: ' num2str(t(4)) ' sec'])
end
disp('-------------------E N D----------------------------')


end %function FRPdriver



function [wintervals]=pass_check_Y(SERflag,A,B,C,D,colinterch)

%  Input: state-space model with diagonal A.  (complex conjugate formulation)
%  
%  Output: matrix with 2 rows. 
%          Each column defines the start and end frequency of non-passive band: [w_start w_end].'.
%          Empty ([]) if there are no violations
%          A violation extending to infinity is denoted by s_end=j*1e16 
%  
%  History: Jan 24, 2006.  Bjorn Gustavsen
% 



%TOL= 1e-6;  %Criterion for deciding whwther an eigenvalue is purely imaginary or not 
wintervals=[];

if SERflag==1 %Must convert from pole-residue to state-space
  Nc=length(D);
  N=length(A);
  tell=0;
  CC=zeros(Nc,Nc*N);
  AA=[]; BB=[]; B=ones(N,1);
  for col=1:Nc
    AA=blkdiag(AA,diag(sparse(A)));
    BB=blkdiag(BB,B);  
    for row=col:Nc
      %tell=tell+1;
      %DD(row,col)=D(tell);
      CC(row,(col-1)*N+1:col*N)=C(row,col,:); 
      CC(col,(row-1)*N+1:row*N)=C(row,col,:);     
    end
  end
  %DD=DD+(DD-diag(diag(DD))).';  
  A=AA; B=BB; C=CC; %D=DD;
end    

Acmplx=A;Bcmplx=B;Ccmplx=C; Dcmplx=D;

if sum(sum(A-diag(diag(A))))==0 %Convert to real-only
  N=length(A);  
  cindex=zeros(1,N);
  for m=1:N 
    if imag(A(m,m))~=0  
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
  n=0;
  for m=1:N
    n=n+1;
    if cindex(m)==1
      a=A(n,n); a1=real(a); a2=imag(a);
      c=C(:,n); c1=real(c); c2=imag(c);
      b=B(n,:); b1=2*real(b); b2=-2*imag(b);
      Ablock=[a1 a2;-a2 a1];
      A(n:n+1,n:n+1)=Ablock;
      C(:,n)=c1;
      C(:,n+1)=c2;
      B(n,:)=b1;
      B(n+1,:)=b2;
    end
  end
end



N=length(A); Nc=length(D);
tell=0;
%B=[];
%for col=1:Nc
%  B=blkdiag(B,ones(N/Nc,1));  
%end
%size(A),size(B),size(C),size(D)

%disp('Checking passivity...')

% %Calculating Hamiltonian matrix:
% Ham11 = A - B*(D+D.')^(-1)*C;
% Ham12 = B*(D+D.')^(-1)*B.';
% Ham21 = -C.'*(D+D.')^(-1)*C;
% Ham22 = -A.' + C.'*(D+D.')^(-1)*B.';
% Ham=[Ham11 Ham12; Ham21 Ham22];

%############
%A=Acmplx; B=Bcmplx;C=Ccmplx;
%############

% %tic
% [Eham]=eig(Ham);
% 
% %[Eham]=eig(Ham,'nobalance');
% %toc
% Eham=flipud(Eham);
% 
% 
% Nc=length(D);
E=zeros(Nc,Nc); %Dummy E-matrix
% 
% %######################
% %Eham=abs(real(Eham))+j*abs(imag(Eham));
% %######################
% 
% Eham_pos=Eham(imag(Eham)>0); %Evals with positive imag. part 
% Eham_pos=Eham_pos+1e-16; %Avoid warning message about division by zero
% ratio=abs(imag(Eham_pos)./real(Eham_pos));
% if spy==1
%   figure(9),
%   h1=semilogy(1./sort(ratio),'b.'); hold on, ylim([1e-16 1e0]);grid on, xlabel('Eigenvalue');ylabel('|Real part./Imag. part|');
%   %h2=semilogy([1 length(ratio)],[TOL,TOL],'r-');
%   title('Hamiltonian matrix eigenvalues')
%   %legend([h1 h2],'Eigenvalue','TOL Hamiltonian',3);
%   drawnow, pause(0.01);hold off
% end  

%#########################################################
if sum(eig(D)==0)>0 %singular D;  
  Ahat=A\eye(N);
  Bhat=-Ahat*B;
  Chat=C*Ahat;
  Dhat=D-C*Ahat*B;
  A=Ahat;B=Bhat;C=Chat;D=Dhat;
end
%!!S1=A*(B/D*C-A);
S1=A*(B*D^(-1)*C-A);
%S1=(B/D*C-A)*A;
%[sqrt(eig(S1)) sqrt(eig(S2))]/(2*pi)

%A=Aold;B=Bold;C=Cold;D=Dold;  
wS1=eig(S1);
wS1=sqrt(wS1);
%wS1/(2*pi)
if sum(eig(Dcmplx)==0)>0 
  wS1=1./wS1;
end  
ind=find(imag(wS1)==0);
wS1=wS1(ind);
sing_w=sort(wS1);


% figure(101);
%         %s_pass=2*pi*i*logspace(0,log10(20),401);I=ones(length(A),1);
%         s_pass=2*pi*i*linspace(0,5,401);I=ones(length(A),1);
%         for k=1:length(s_pass)
%           Y   =Ccmplx*diag( ( s_pass(k)*I -diag(Acmplx) ).^(-1) )*Bcmplx +Dcmplx +s_pass(k)*E;
%           EE1(:,k)=eig(real(Y)); 
%         end 
%         plot(s_pass/(2*pi*i),EE1);grid on
%         save koko A B C D

if length(sing_w)==0, sing_w=[]; intervals=[]; return, end

A=Acmplx; B=Bcmplx;C=Ccmplx; D=Dcmplx;
%fS1=wS1/(2*pi)
%#########################################################

%#####################
 %[dum1,ind]=sort(ratio);
 %dum2=Eham_pos(ind);
 %[dum2 dum1]
%#####################


% singulars=Eham_pos(ratio>1/TOL);
% %singulars=Eham_pos(1./ratio<TOL);
% 
% if length(singulars)==0, sing_w=[]; intervals=[]; return, end
% sing_w=imag(singulars); %Crossing frequencies [rad/sec]
% % sing_w=sort(sing_w);
  
%Establising frequency list at midpoint of all bands defined by sing_w:
midw=zeros(1+length(sing_w),1); 
midw(1)=sing_w(1)/2;
midw(end)=2*sing_w(end);
for k=1:length(sing_w)-1
  midw(k+1)=(sing_w(k)+sing_w(k+1))/2;
end  

%Checking passivity at all midpoints:
for k=1:length(midw);
  sk=i*midw(k);
  G=real(fitcalcABCDE(sk,diag(A),B,C,D,E));
  EE(:,k)=eig(G); 
  if any(EE(:,k)<0)
    viol(k)=1;
  else
    viol(k)=0;
  end  
end


%Establishing intervals for passivity violations: 
intervals=[];
for k=1:length(midw)
  if viol(k)==1
    if k==1
      intervals=[intervals [0 sing_w(1)].'];  %The first violations starts at DC
    elseif k==length(midw)
      intervals=[intervals [sing_w(k-1)  1e16].'];  %The last violation extends to infinite frequency      
    else
      intervals=[intervals [sing_w(k-1) sing_w(k)].']; 
    end  
  end
end



if isempty(intervals), wintervals=intervals;return,end
%Collapsing overlapping bands:
tell=0;
killindex=0;
for k=2:1:length(intervals(1,:))
  if intervals(2,k-1)==intervals(1,k) %An overlap exists
    tell=tell+1;
    intervals(2,k-1)=intervals(2,k); %Extending interval
    intervals(:,k)=intervals(:,k-1); %Copying interval 
    killindex(tell)=k-1;
  end
end

if killindex~=0
  intervals(:,killindex)=[];
end  
wintervals=intervals;

%disp('Done')




end %function pass_check_Y




%S-PARAMETER CASE:
function [wintervals]=pass_check_S(SERflag,A,B,C,D,TOL,spy,colinterch)

%  Input: state-space model with diagonal A.  (complex conjugate formulation)
%  
%  Output: matrix with 2 rows. 
%          Each column defines the start and end frequency of non-passive band: [w_start w_end].'.
%          Empty ([]) if there are no violations
%          A violation extending to infinity is denoted by s_end=j*1e16 
%  
%  History: Jan 24, 2006.  Bjorn Gustavsen
% 



%TOL= 1e-6;  %Criterion for deciding whwther an eigenvalue is purely imaginary or not 
wintervals=[];

if SERflag==1 %Must convert from pole-residue to state-space
  Nc=length(D);
  N=length(A);
  tell=0;
  CC=zeros(Nc,Nc*N);
  AA=[]; BB=[]; B=ones(N,1);
  for col=1:Nc
    AA=blkdiag(AA,diag(sparse(A)));
    BB=blkdiag(BB,B);  
    for row=col:Nc
      %tell=tell+1;
      %DD(row,col)=D(tell);
      CC(row,(col-1)*N+1:col*N)=C(row,col,:); 
      CC(col,(row-1)*N+1:row*N)=C(row,col,:);     
    end
  end
  %DD=DD+(DD-diag(diag(DD))).';  
  A=AA; B=BB; C=CC; %D=DD;
end    

Acmplx=A;Bcmplx=B;Ccmplx=C;

if sum(sum(A-diag(diag(A))))==0 %Convert to real-only
  N=length(A);  
  cindex=zeros(1,N);
  for m=1:N 
    if imag(A(m,m))~=0  
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
  n=0;
  for m=1:N
    n=n+1;
    if cindex(m)==1
      a=A(n,n); a1=real(a); a2=imag(a);
      c=C(:,n); c1=real(c); c2=imag(c);
      b=B(n,:); b1=2*real(b); b2=-2*imag(b);
      Ablock=[a1 a2;-a2 a1];
      A(n:n+1,n:n+1)=Ablock;
      C(:,n)=c1;
      C(:,n+1)=c2;
      B(n,:)=b1;
      B(n+1,:)=b2;
    end
  end
end




N=length(A); Nc=length(D);
tell=0;
%B=[];
%for col=1:Nc
%  B=blkdiag(B,ones(N/Nc,1));  
%end
%size(A),size(B),size(C),size(D)

%disp('Checking passivity...')

%Calculating Hamiltonian matrix:
I=eye(length(D));R=D.'*D-I; 
invR=R^(-1);
% S=R;
% invS=S^(-1);
% Ham11 = A - B*invR*D.'*C;
% Ham12 = -B*invR*B.';
% Ham21 = C.'*invR*C;
% Ham22 = -A.' + C.'*D*invR*B.';
% Ham=[Ham11 Ham12; Ham21 Ham22];
% %############
% A=Acmplx; B=Bcmplx;C=Ccmplx;
% %############
% [Eham]=eig(Ham);
% Eham=flipud(Eham);
% Nc=length(D);
% E=zeros(Nc,Nc); %Dummy E-matrix
% Eham_pos=Eham(imag(Eham)>0); %Evals with positive imag. part 
% Eham_pos=Eham_pos+1e-16; %Avoid warning message about division by zero
% ratio=abs(imag(Eham_pos)./real(Eham_pos));
% singulars=Eham_pos(ratio>1/TOL)


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%E=A - B*invR*D.'*C;
%F=-B*invR*C;
%P=(E+F)*(E-F);
P=(A-B*(D-I)^(-1)*C)*(A-B*(D+I)^(-1)*C);
wS1=sqrt(eig(P));
ind=find(real(wS1)==0);
wS1=imag(wS1(ind));
singulars=j*sort(wS1);
A=Acmplx; B=Bcmplx;C=Ccmplx;
E=zeros(Nc,Nc); %Dummy E-matrix
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


if length(singulars)==0, sing_w=[]; intervals=[]; return, end
sing_w=imag(singulars); %Crossing frequencies [rad/sec]
sing_w=sort(sing_w);
   
%Establising frequency list at midpoint of all bands defeined by sing_w:
midw=zeros(1+length(sing_w),1); 
midw(1)=sing_w(1)/2;
midw(end)=2*sing_w(end);
for k=1:length(sing_w)-1
  midw(k+1)=(sing_w(k)+sing_w(k+1))/2;
end  

%Checking passivity at all midpoints:
oldU=[];
for k=1:length(midw);
  sk=i*midw(k);
  Y=(fitcalcABCDE(sk,diag(A),B,C,D,E));
  if colinterch==0  
    EE(:,k)=svd(Y,0); 
  else 
    [U,S,V]=svd(Y,0);   
    [U,S,V]=interchsvd(U,oldU,S,V,Nc,k); oldU=U;   
    EE(:,k)=diag(S);   
  end 
  if any(EE(:,k)>1)
    viol(k)=1;
  else
    viol(k)=0;
  end  
end


%Establishing intervals for passivity violations: 
intervals=[];
for k=1:length(midw)
  if viol(k)==1
    if k==1
      intervals=[intervals [0 sing_w(1)].'];  %The first violations starts at DC
    elseif k==length(midw)
      intervals=[intervals [sing_w(k-1)  1e16].'];  %The last violation extends to infinite frequency      
    else
      intervals=[intervals [sing_w(k-1) sing_w(k)].']; 
    end  
  end
end

if isempty(intervals), wintervals=intervals;return,end
%Collapsing overlapping bands:
tell=0;
killindex=0;
for k=2:1:length(intervals(1,:))
  if intervals(2,k-1)==intervals(1,k) %An overlap exists
    tell=tell+1;
    intervals(2,k-1)=intervals(2,k); %Extending interval
    intervals(:,k)=intervals(:,k-1); %Copying interval 
    killindex(tell)=k-1;
  end
end

if killindex~=0
  intervals(:,killindex)=[];
end  
wintervals=intervals;

%disp('Done')


end %function pass_check_ham














 

 

