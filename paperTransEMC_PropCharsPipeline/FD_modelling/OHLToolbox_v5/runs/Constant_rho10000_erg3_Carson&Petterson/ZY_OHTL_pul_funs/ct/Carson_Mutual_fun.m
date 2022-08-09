function [Z12_Carson] = Carson_Mutual_fun(num,h,d,con,e_g,omega,rad_ex)
%
% Carson_Mutual computes the earth-return mutual impedance of an overhead
% conductor above earth. Since permittivity is included it can evaluate 
% also Sunde's result. One of the parameters can also be vector of values.
%
% Author : T.P. Theodoulidis
% Date   : 21 March 2015
%
% Arguments
% h1 : conductor-1 height [m]
% h2 : conductor-2 height [m]
% x  : distance between conductors [m]
% f  : frequency [Hz]
% con: earth conductivity [S/m]
% er : relative dielectric permittivity (er=0 gives Carson's result)
% Z12: earth-return mutual-impedance [Ohm]
% CI : Carson's integral
%
% External routines called : StruveH1Y1 available from Matlab Central
%

%% create empty matrices for Z12 and CI


Z12_Carson=zeros(num,num);
CI_Carson=zeros(num,num);

%%constants

m0=pi*4e-7;
%e0=8.8541878176e-12;

%% filling the matrices 
for x=1:1:num
    for y=1:1:num
        
%             omega=2*pi*f;

         if x~=y
             %mutual impedances
             e_g=0;
            k=sqrt(1i*omega*m0.*(con+1i*omega*e_g));
            H=h(1,x)+h(1,y);
            d1=sqrt((h(1,x)-h(1,y))^2+d(x,y)^2);
            d2=sqrt(H^2+d(x,y)^2);
            u1=k*(H-1i*d(x,y));
            u2=k*(H+1i*d(x,y));
            %
            CI_Carson(x,y) = (pi/2./u1).*(StruveH1Y1(u1)-2/pi./u1)+...
                 (pi/2./u2).*(StruveH1Y1(u2)-2/pi./u2);
            %
            Z12_Carson(x,y) = (1i.*omega*m0/(2*pi))*(log(d2/d1)+CI_Carson(x,y));
%          else
             
             
%              % self impedances
%                 k=sqrt(1i*omega*m0.*(con+1i*omega*e_g));
%                 d1=rad_ex(x,x);
%                 d2=2*h(1,x);
%                 u=k*(2*h(1,x));
%                 %
%                 CI_Carson(x,y)= 2*(pi/2./u).*(StruveH1Y1(u)-2/pi./u);
%                 %
%                 Z12_Carson(x,y) = 1i.*omega*m0/(2*pi)*(log(d2/d1)+CI_Carson(x,y));

%                  
         end

    end
end

