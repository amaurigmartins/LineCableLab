function []=netgen_ATP(SER,NOD,fname);

% []=netgen_ATP(SER,NOD,fname);
%
% This file is part of the Matrix Fitting Toolbox, v1.
% Filename: netgen_ATP.m
% Package: Matrix_Fitting_Toolbox_1.zip.
% Programmed by B. Gustavsen. October 08, 2008.
%
% PURPOSE:
% The program generates electrical network for a matrix Y whose 
% elements have been fitted with an identical set of poles (pole-residue model).
% The electrical network is written to file in a format which can be read by ATP.
%
% INPUT:
% SER: structure containing pole residue model (produced by poleresiduefit.m)
%      SER.R     (Nc,Nc,N) :Residue matrix
%      SER.poles (1,N) 
%      SER.D     (Nc,Nc)
%      SER.E     (Nc,Nc)
%
% NOD: Single character name that defines the node name in the ATP-data file.
%      Valid values: 'A','B',...'Z'  
%   
%      For instance, choosing 'X' results in that the model is connected to node names  
%      X____1, X____2, X____3,...
%
% fname: string defining the file name 



disp('****Creating equivalent circuit for ATP...')


SERA=SER.poles;
SERC=SER.R;
SERD=SER.D;
SERE=SER.E;


Nc=length(SERD);
N =length(SERA);

%Avoiding problems with zero D and E,
SERD=SERD+1e-16 +1e-16*eye(Nc);
SERE=SERE+1e-16 +1e-16*eye(Nc);
for m=1:length(SERC)
  SERC(:,:,m)=SERC(:,:,m)+1e-16 +1e-16*eye(Nc);
end  


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

%=================================

%fid1=fopen('RLC.out','w');
fid1=fopen(fname,'wt');
fprintf(fid1,'$VINTAGE,1\n');


fprintf(fid1,'C <BUS1><BUS2><BUS3><BUS4><   OHM        ><   milliH     ><   microF     >\n');
fprintf(fid1,'C \n');

fprintf(fid1,'$VINTAGE,1\n');
  
fprintf(fid1,'C <BUS1><BUS2><BUS3><BUS4><   OHM        ><   milliH     ><   microF     >\n');
fprintf(fid1,'C \n');

for row=1:Nc
 for col=row:Nc

  fprintf(fid1,'C (%1.0f,%1.0f)\n',[row col]);

  if row==col %diagonal element
    dum=SERD(row,:);dum(col)=0;
    D=SERD(row,col)+sum(dum);
      R0=1/D;
    %!!fprintf(fid1,'  A____');fprintf(fid1,'%1.0f',row');
      fprintf(fid1,'%s%s%s','  ',NOD,'____');fprintf(fid1,'%1.0f',row');           
      fprintf(fid1,'                  ');
      fprintf(fid1,'%16.8e\n',R0);
    dum=SERE(row,:);dum(col)=0;
    C0=SERE(row,col)+sum(dum);
    %!!fprintf(fid1,'  A____');fprintf(fid1,'%1.0f',row');
      fprintf(fid1,'%s%s%s','  ',NOD,'____');fprintf(fid1,'%1.0f',row');            
      fprintf(fid1,'                                                  ');
      fprintf(fid1,'%16.8e\n',1e6*C0);

    for m=1:N
      M=m+(m<0)*2^32; %Convert to hex
      if cindex(m)==0 %real pole
        a1=SERA(m);
        dum=squeeze(SERC(row,:,m));dum(col)=0;
        c1=SERC(row,col,m)+sum(dum);
        L1=1/c1;
        R1=-a1/c1;

       %!!fprintf(fid1,'  A____');fprintf(fid1,'%1.0f',row');
          fprintf(fid1,'%s%s%s','  ',NOD,'____');fprintf(fid1,'%1.0f',row');               
          fprintf(fid1,'                  ');
          fprintf(fid1,'%16.8e',R1);
          fprintf(fid1,'%16.8e\n',1000*L1);
      elseif cindex(m)==1 %complex pole, 1st part
        a1=real(squeeze(SERA(m)));
        a2=imag(squeeze(SERA(m)));  
        dum=squeeze(SERC(row,:,m));dum(col)=0;
        dum=SERC(row,col,m)+sum(dum);
        c1=real(dum);
        c2=imag(dum);
        L=1/(2*c1);
        dum=c1*a1+c2*a2;
        R=(-2*a1+2*dum*L)*L;
        C=(a1^2+a2^2+2*dum*R)*L; C=1/C;
        G=-2*dum*C*L;
       %!! fprintf(fid1,'  A____');fprintf(fid1,'%1.0f',row');
          fprintf(fid1,'%s%s%s','  ',NOD,'____');fprintf(fid1,'%1.0f',row');        
       %!!fprintf(fid1,'A'); fprintf(fid1,'%2X',M);fprintf(fid1,'__');fprintf(fid1,'%1.0f',row');
          fprintf(fid1,'%s%2X',NOD,M);fprintf(fid1,'__');fprintf(fid1,'%1.0f',row');   

          fprintf(fid1,'            ');
          fprintf(fid1,'%16.8e',R);
          fprintf(fid1,'%16.8e\n',1000*L);

       %!!fprintf(fid1,'  A'); fprintf(fid1,'%2X',M);fprintf(fid1,'__');fprintf(fid1,'%1.0f',row');
          fprintf(fid1,'%s%s%2X','  ',NOD,M);fprintf(fid1,'__');fprintf(fid1,'%1.0f',row');   

          fprintf(fid1,'                  ');
          fprintf(fid1,'%16.8e\n',1/G);
       %!!fprintf(fid1,'  A'); fprintf(fid1,'%2X',M);fprintf(fid1,'__');fprintf(fid1,'%1.0f',row');
          fprintf(fid1,'%s%s%2X','  ',NOD,M);fprintf(fid1,'__');fprintf(fid1,'%1.0f',row');   
          fprintf(fid1,'                                                  ');
          fprintf(fid1,'%16.8e\n',1e6*C);    

      end %if cindex(m)==

     end %for m=1:N

   else %row~=col (off-diagonal element)

    C0=-SERE(row,col);
      R0=-1/SERD(row,col); 
   %!!fprintf(fid1,'  A____');fprintf(fid1,'%1.0f',row');
   %!!fprintf(fid1,'A____');fprintf(fid1,'%1.0f',col');
      fprintf(fid1,'%s%s%s','  ',NOD,'____');fprintf(fid1,'%1.0f',row');      
      fprintf(fid1,'%s%s',NOD,'____');fprintf(fid1,'%1.0f',col');         
      fprintf(fid1,'           ');
      fprintf(fid1,' %16.8e\n',R0);
   %!!fprintf(fid1,'  A____');fprintf(fid1,'%1.0f',row');
   %!!fprintf(fid1,'A____');fprintf(fid1,'%1.0f',col');
      fprintf(fid1,'%s%s%s','  ',NOD,'____');fprintf(fid1,'%1.0f',row');      
      fprintf(fid1,'%s%s',NOD,'____');fprintf(fid1,'%1.0f',col');           
      fprintf(fid1,'                                            ');
      fprintf(fid1,'%16.8e\n',1e6*C0);

    for m=1:N
      M=m+(m<0)*2^32; %Convert to hex
      if cindex(m)==0 %real pole
        a1=SERA(m);
        c1=-SERC(row,col,m);
        L1=1/c1;
        R1=-a1/c1;
       %!!fprintf(fid1,'  A____');fprintf(fid1,'%1.0f',row');
       %!!fprintf(fid1,'A____');fprintf(fid1,'%1.0f',col');
          fprintf(fid1,'%s%s%s','  ',NOD,'____');fprintf(fid1,'%1.0f',row');      
          fprintf(fid1,'%s%s',NOD,'____');fprintf(fid1,'%1.0f',col');    
          fprintf(fid1,'            ');
          fprintf(fid1,'%16.8e',R1);
          fprintf(fid1,'%16.8e\n',1000*L1);
      elseif cindex(m)==1 %complex pole, 1st part

        a1=real(squeeze(SERA(m)));
        a2=imag(squeeze(SERA(m)));  
        c1=-real(SERC(row,col,m));
        c2=-imag(SERC(row,col,m));
        L=2*c1; L=1/L;
        dum=c1*a1+c2*a2;
        R=(-2*a1+2*dum*L)*L;
        C=(a1^2+a2^2+2*dum*R)*L; C=1/C;
        G=-2*dum*C*L;
       %!!fprintf(fid1,'  A____');fprintf(fid1,'%1.0f',row');
          fprintf(fid1,'%s%s%s','  ',NOD,'____');fprintf(fid1,'%1.0f',row');          
       %!!fprintf(fid1,'A'); fprintf(fid1,'%2X',M);fprintf(fid1,'_');fprintf(fid1,'%1.0f',row,col');
          fprintf(fid1,'%s%2X',NOD,M);fprintf(fid1,'_');fprintf(fid1,'%1.0f',row,col');   

          fprintf(fid1,'            '); 
          fprintf(fid1,'%16.8e',R);
          fprintf(fid1,'%16.8e\n',1000*L);

       %!!fprintf(fid1,'  A'); fprintf(fid1,'%2X',M);fprintf(fid1,'_');fprintf(fid1,'%1.0f',row,col');
          fprintf(fid1,'%s%s%2X','  ',NOD,M);fprintf(fid1,'_');fprintf(fid1,'%1.0f',row,col');   

      %!! fprintf(fid1,'A____');fprintf(fid1,'%1.0f',col');
          fprintf(fid1,'%s%s',NOD,'____');fprintf(fid1,'%1.0f',col');   

          fprintf(fid1,'            ');
          fprintf(fid1,'%16.8e\n',1/G);
      %!! fprintf(fid1,'  A'); fprintf(fid1,'%2X',M);fprintf(fid1,'_');fprintf(fid1,'%1.0f',row,col');
          fprintf(fid1,'%s%s%2X','  ',NOD,M);fprintf(fid1,'_');fprintf(fid1,'%1.0f',row,col');   

      %!! fprintf(fid1,'A____');fprintf(fid1,'%1.0f',col');
          fprintf(fid1,'%s%s',NOD,'____');fprintf(fid1,'%1.0f',col');              
          fprintf(fid1,'                                            ');
          fprintf(fid1,'%16.8e\n',1e6*C);

      end %if cindex(m)==

     end %for m=1:N

   end 

 end %for col=
end %for row=


fprintf(fid1,'$VINTAGE,0\n');
fclose(fid1);





 