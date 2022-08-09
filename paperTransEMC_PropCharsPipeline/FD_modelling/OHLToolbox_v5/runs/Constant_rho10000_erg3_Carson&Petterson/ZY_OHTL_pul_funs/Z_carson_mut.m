function [Z12] = Z_carson_mut(h,d,omega,con,er,ord)
%Z_carson_mut(h1,h2,x,f,con,er)
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
m0=pi*4e-7;
e0=8.8541878176e-12;
%omega=2*pi*f;

Z12=zeros(ord,ord);
CI=zeros(ord,ord);

for x=1:1:ord
    for y=1:1:ord
       
       if x~=y 
        k=sqrt(1i*omega*m0.*(con+1i*omega*e0*er));
        H=h(1,x)+h(1,y);
        d1=sqrt((h(1,x)-h(1,y))^2+d(x,y)^2);
        d2=sqrt(H^2+d(x,y)^2);
        u1=k*(H-1i*d(x,y));
        u2=k*(H+1i*d(x,y));
        %
        CI(x,y) = (pi/2./u1).*(StruveH1Y1(u1)-2/pi./u1)+...
             (pi/2./u2).*(StruveH1Y1(u2)-2/pi./u2);
        %
        Z12(x,y) = 1i.*omega*m0/(2*pi)*(log(d2/d1)+CI(x,y));
       end
    end
end
%
