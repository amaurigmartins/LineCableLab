%function [Z1,CI] = Carson_Self(h,r,f,con,er)
function [Z12_Carson] = Carson_Self_fun(num,h,d,con,e_g,omega,rad_ex)
%
% Carson_Self computes the earth-return self impedance of an overhead
% conductor above earth. Since permittivity is included it can evaluate 
% also Sunde's result.  One of the parameters can also be vector of values.
%
% Author : T.P. Theodoulidis
% Date   : 21 March 2015
%
% Arguments
% h  : conductor height [m]
% r  : conductor radius [m]
% f  : frequency [Hz]
% con: earth conductivity [S/m]
% er : relative dielectric permittivity (er=0 gives Carson's result)
% Z1 : earth-return self-impedance [Ohm]
% CI : Carson's integral
%
% External routines called : StruveH1Y1 available from Matlab Central
%
m0=pi*4e-7;
Z12_Carson=zeros(num,num);
CI_Carson=zeros(num,num);


%% filling the matrices 
for x=1:1:num
   
        e_g=0;
        k=sqrt(1i*omega*m0.*(con+1i*omega*e_g));
        d1=rad_ex(x,x);
        d2=2*h(1,x);
        u=k*(2*h(1,x));
        %
        CI_Carson(x,x) = 2*(pi/2./u).*(StruveH1Y1(u)-2/pi./u);
        %
        Z12_Carson(x,x) = 1i.*omega*m0/(2*pi)*(log(d2/d1)+CI_Carson(x,x));
    
end
%
