function [length,Ncon,soil,h,d,Geom]=LineData_fun()
% Line geometry
% 1 column -- number of phase (set to 0 for Kron reduction)
% 2 column -- x position of each conduntor in meters
% 3 column -- y position of each coductor in meters
% 4 column -- internal radii of each conductor
% 5 column -- external radii of each conductor
% 6 column -- resistivity of the conductor
% 7 column -- permeability of the conductor
% 8 column -- external radii of insulation
% 9 column -- relative permeability of insulation
% 10 column -- relative permittivity of insulation
% 11 column -- line length in m
length = 1000.000000;
Geom=[ ...
     1  -6.600000  13.500000  0.004630  0.012570  7.1221e-08  1.000000     NaN     NaN     NaN  length
     2  0.000000  13.500000  0.004630  0.012570  7.1221e-08  1.000000     NaN     NaN     NaN  length
     3  6.600000  13.500000  0.004630  0.012570  7.1221e-08  1.000000     NaN     NaN     NaN  length
     4  10.000000  1.000000  0.124500  0.127000  2.844e-07  250.000000  0.227000  1.000000  3.000000  length
     5  -4.650000  17.600000  0.000000  0.004765  2.46925e-07  1.000000     NaN     NaN     NaN  length
     6  4.650000  17.600000  0.000000  0.004765  2.46925e-07  1.000000     NaN     NaN     NaN  length
];
Ncon = size(Geom,1);
% Variables
e0=8.854187817e-12;
m0=4*pi*1e-7;
% Rearrange data matrix to perform bundle reduction
Geom = reorderGeoMatrix(Geom);
% Height of line calculation
[h]=height_fun(Ncon,Geom);
% Distances calculations
[d]=distance_fun(Ncon,Geom);
% Earth electric parameters
soil.erg = e0*3.000000;
soil.m_g = m0*1.000000;
soil.sigma_g = 1/1000.000000;
end