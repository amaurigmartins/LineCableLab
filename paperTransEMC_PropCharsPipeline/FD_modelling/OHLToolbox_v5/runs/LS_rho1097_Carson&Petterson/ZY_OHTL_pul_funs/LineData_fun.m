function [length,Ncon,soil,h,d,Geom]=LineData_fun()
% Line Geometry
%
% 1 column -- number of conductor
% 2 column -- x position of each conduntor in meters
% 3 column -- y position of each coductor in meters
% 4 column -- internal radii of each conductor
% 5 column -- external radii of each conductor
% 6 column -- resistivity of the aluminum
% 7 column -- permeability of the conductor
% 8 column -- external radii of insulation
% 9 column -- relative permeability of insulation
% 10 column -- relative permittivity of insulation
% 11 column -- line length in m


% Geom = [1   0.0     20   0.00463  0.01257  7.1221e-8   1     nan    nan   nan   1
%         2   10.0     20   0.00463  0.01257  7.1221e-8   1     nan    nan   nan   1];
    
Geom = [1  -6.6     13.5   0.00463  0.01257  7.1221e-8   1     nan    nan   nan   1
        2   0.0     13.5   0.00463  0.01257  7.1221e-8   1     nan    nan   nan   1
        3   6.6     13.5   0.00463  0.01257  7.1221e-8   1     nan    nan   nan   1
        4  10.0     01.0   0.12450  0.12700  2.8444e-7   250   0.227    1     3   1       % P46
        5  -4.65    17.6   0.00000  0.004765 2.46925E-7   1     nan    nan   nan  1
        6   4.65    17.6   0.00000  0.004765 2.46925E-7   1     nan    nan   nan   1];     % P66       
%        4  10.0     01.0   0.12450  0.12700  2.8444e-7   250   nan    nan   nan   1];    % P44       
%        4  10.0     01.0   0.12450  0.12700  2.8444e-7   1   nan    nan   nan   1];      % P43       
%        4  10.0     01.0   0.00463  0.01257  7.1221e-8   1   nan    nan   nan   1];      % P42       
%        4  10.0     13.5   0.00463  0.01257  7.1221e-8   1   nan    nan   nan   1];      % P41           
length  = Geom(1,11);                                     % Line length
Ncon    = Geom(max(Geom(:,1)),1);                        % Number of conductors

% Variables
%e0=8.854187817e-12;  % Farads/meters
m0=4*pi*1e-7;        % Henry's/meters


%Height of Line Calculation
[h]=height_fun(Ncon,Geom);

%Distance between conductor calculation
[d]=distance_fun(Ncon,Geom);

% Earth Electric Parameters
soil.erg=1; %relative permittivity of earth
mrg=1;          %relative permeability of earth
soil.m_g=m0*mrg;
rho=1097;%permeability of earth
soil.sigma_g=1/rho; %conductivity of earth

