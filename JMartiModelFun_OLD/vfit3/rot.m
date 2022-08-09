function [OUT]=rot(S);    %03.04.2001, Bjorn G.

% [OUT]=rot(S);
%
% This file is part of the Matrix Fitting Toolbox, v1.
% Filename: rot.m
% Package: Matrix_Fitting_Toolbox_1.zip.
% Programmed by B. Gustavsen. October 08, 2008.


Nc=length(S);
SA=zeros(Nc); SB=SA; scale1=zeros(1,Nc); scale2=scale1; err1=scale1;err2=scale1;
% E-vector #col :      
      for col=1:Nc
% Calculate min/max-value for square sum (error) of imag.parts :
          numerator(col)=0.0;
          denominator(col)=0.0;
        for j=1:Nc
          numerator(col)=numerator(col)+imag(S(j,col))*real(S(j,col));           %   Q2(j,i)*Q1(j,i)
          denominator(col)=denominator(col)+(real(S(j,col)))^2 - (imag(S(j,col)))^2; %   Q1(j,i)**2 - Q2(j,i)**2
        end

        numerator(col)=-2.0*numerator(col);
        ang(col)=0.5*atan2(numerator(col),denominator(col));

        scale1(col)=cos(ang(col)) +i*sin(ang(col));
        %11.07.2001 B.G. 
 %       scale2(col)=scale1(col) + pi/2.0;   
        scale2(col)=cos(ang(col)+pi/2) +i*sin(ang(col)+pi/2);
% Deciding which solution (1,2) will produce the smallest error :

        for j=1:Nc
          SA(j,col)=S(j,col)*scale1(col);
          SB(j,col)=S(j,col)*scale2(col);
        end

% Square sum (error) of solution :
           aaa=0.0;
           bbb=0.0;    
           ccc=0.0;    
        for j=1:Nc
          aaa = aaa +  (imag(SA(j,col)))^2;                % Q2A(j,i)**2
          bbb = bbb +  real(SA(j,col))*imag(SA(j,col));      % Q1A(j,i)*Q2A(j,i)
          ccc = ccc +  (real(SA(j,col)))^2;                % Q1A(j,i)**2
        end
          err1(col)=aaa*cos(ang(col))^2 + bbb*sin(2.0d0*ang(col))+ ccc*sin(ang(col))^2;

% Square sum (error) of solution #2 :
           aaa=0.0;    
           bbb=0.0;    
           ccc=0.0;    
        for j=1:Nc
          aaa = aaa +  (imag(SB(j,col)))^2;                % Q2A(j,i)**2
          bbb = bbb +  real(SB(j,col))*imag(SB(j,col));      % Q1A(j,i)*Q2A(j,i)
          ccc = ccc +  (real(SB(j,col)))^2;                % Q1A(j,i)**2
        end
          err2(col)=aaa*cos(ang(col))^2 + bbb*sin(2.0d0*ang(col))+ ccc*sin(ang(col))^2;%+ ccc*sin(ang(col))^2;
                                  
% Picking the solution (1,2) with the smallest square sum :
        if(err1(col)<err2(col)) 
          scale(col)=scale1(col);
        else
          scale(col)=scale2(col);
        end

% Rotating e-vector #col :
        S(:,col)=S(:,col).*scale(col);
end

      OUT=S;