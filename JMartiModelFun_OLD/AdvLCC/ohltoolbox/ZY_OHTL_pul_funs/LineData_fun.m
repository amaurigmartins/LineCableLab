function [length,Ncon,r,rad_ex,rad_in,erg,mrg,sigma_g,sigma_w,e0,m0,m_g,h,d]=LineData_fun()
% Line Geometry
%
% 1 column -- number of conductor
% 2 column -- x position of each conduntor in meters
% 3 column -- y position of each coductor in meters
% 4 column -- radii of each conductor
% 5 column -- resistivity of the aluminum
% 6 column -- relative permitivity
% 7 column -- line length in m

Geom = [1   0     8   0.00318  2.826e-8   10   100
        2   1.5   8   0.00318  2.826e-8   10   100
        3   2.3   8   0.00318  2.826e-8   10   100];
             
length  = Geom(1,7);                                     % Line length
Ncon    = Geom(max(Geom(:,1)),1);                        % Number of conductors

r=0.00862;      % external radius of conductor
rad_ex=r;
rad_in=Geom(1,4); % internal radius of conductor

% Electric Parameters
erg=Geom(1,6);         %relative permittivity of earth
mrg=1;          %relative permeability of earth
sigma_g=1/1000; %conductivity of earth
sigma_w=1/3.86063e-008; %conductivity of conductor

% Variables
e0=8.854187817e-12;  % Farads/meters
m0=4*pi*1e-7;        % Henry's/meters
m_g=m0*mrg;

%Height of Line Calculation
for i=1:1:Ncon
    [h]=height_fun(Ncon,i,Geom);
end

%Distance between conductor calculation
[d]=distance_fun(Ncon,Geom);
