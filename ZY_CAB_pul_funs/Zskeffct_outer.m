% Zskeffct_outer.m
function [zin]=Zskeffct_outer(radius_ex,radius_in,sigma_c,mu_r,omega,varargin)
% Calculation of outer surface impedance of a tubular conductor

% Parse optional input
p = inputParser;
addParameter(p, 'SimplifiedFormula', false);
parse(p, varargin{:});
SimplifiedFormula = p.Results.SimplifiedFormula;

% Constant Variables
m0=4*pi*1e-7;
mu_c=m0*mu_r;
TOL=1e-6;

m=sqrt(1i.*omega.*mu_c.*sigma_c); %reciprocal of skin depth

% Approximated skin effect
if (~radius_in)
    radius_in=eps;
end

if SimplifiedFormula
    if radius_in < TOL
        cothTerm = coth(m.*radius_ex.*0.733);
    else
        cothTerm = coth(m.*(radius_ex-radius_in));
    end
    Z1 = (m./sigma_c)./(2.*pi.*radius_ex).*cothTerm;

    if radius_in < TOL
        Z2 = 0.3179./(sigma_c.*pi.*radius_ex.^2);
    else
        Z2 = 1./(sigma_c.*2.*pi.*radius_ex.*(radius_in + radius_ex));
    end
    zin=Z1+Z2;

else
    w_out=m.*radius_ex;
    w_in=m.*radius_in;

    s_in = exp(abs(real(w_in)) - w_out);
    s_out = exp(abs(real(w_out )) - w_in);
    sc = s_in./s_out; %should be applied to all besseli() involving w_in

    N = besseli(0, w_out, 1).*besselk(1, w_in, 1) + sc.*besselk(0, w_out, 1).*besseli(1, w_in, 1);
    D = besseli(1, w_out, 1).*besselk(1, w_in, 1) - sc.*besselk(1, w_out, 1).*besseli(1, w_in, 1);

    % zin = (m./(2.*pi.*radius_ex.*sigma_c)).*(N./D);
    zin = (1i*omega*mu_c/(2*pi)).*(1/w_out).*(N./D);
end

end
