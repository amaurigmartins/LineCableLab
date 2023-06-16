% skeffct_tb_fun.m
function [zin]=skeffct_tb_fun_ct(radius_ex,radius_in,sigma_w,mrw,omega)
% Calculation of Skin Effect for Tubular Conductor

% Constant Variables
m0=4*pi*1e-7;
m_w=m0*mrw;

currmfile = mfilename('fullpath');
currPath = currmfile(1:end-length(mfilename()));

mparam=sqrt(1i.*omega.*m_w.*sigma_w);

% Approximated skin effect
if (~radius_in)
    radius_in=eps;
end

param=mparam.*radius_ex;
param2=mparam.*radius_in;
Zin_tmp1=mparam./(2*pi*radius_ex.*sigma_w);
Zin_tmp2=sym(besseli(0,param).*besselk(1,param2)+besselk(0,param).*besseli(1,param2));
D=sym(besseli(1,param).*besselk(1,param2)-besselk(1,param).*besseli(1,param2));

zin=double(vpa(Zin_tmp1./D.*Zin_tmp2));
global pythoncall

if isnan(zin)
    pycalc_zin=@(Zin_tmp1,param,param2) system([pythoncall ' ' fullfile(currPath,'pycalc_zin.py') ' ' num2str(real(Zin_tmp1)) ' ' num2str(imag(Zin_tmp1)) ' ' num2str(real(param)) ' ' num2str(imag(param)) ' ' num2str(real(param2)) ' ' num2str(imag(param2))]);
    [~,out]=pycalc_zin(Zin_tmp1,param,param2);
    zin=str2num(out);
    if isempty(zin) %this means that python is not installed or not part of the %PATH%
        zin=0;
    end
end

end