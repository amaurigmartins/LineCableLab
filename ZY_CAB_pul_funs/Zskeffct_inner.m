% Zskeffct_inner.m
function [zin]=Zskeffct_inner(radius_ex,radius_in,sigma_c,mu_r,omega, varargin)
% Calculation of inner surface impedance of a tubular conductor

% Parse optional input
p = inputParser;
addParameter(p, 'SimplifiedFormula', false);
parse(p, varargin{:});
SimplifiedFormula = p.Results.SimplifiedFormula;

% Constant Variables
m0=4*pi*1e-7;
mu_c=m0*mu_r;


m=sqrt(1i.*omega.*mu_c.*sigma_c); %reciprocal of skin depth

% Approximated skin effect
if (~radius_in)
    radius_in=eps;
end

if SimplifiedFormula
    Z1 = (m./sigma_c)./(2.*pi.*radius_in).*coth(m.*(radius_ex-radius_in));
    Z2 = 1./(2.*pi.*radius_in.*(radius_in+radius_ex).*sigma_c);
    zin = Z1 + Z2;
else
    w_out=m.*radius_ex;
    w_in=m.*radius_in;

    s_in = exp(abs(real(w_in)) - w_out);
    s_out = exp(abs(real(w_out )) - w_in);
    sc = s_in./s_out; %should be applied to all besseli() involving w_in

    N = sc.*besseli(0, w_in, 1) .*besselk(1, w_out, 1) +     besselk(0, w_in, 1) .*besseli(1, w_out, 1);
    D =     besseli(1, w_out, 1).*besselk(1, w_in, 1)  - sc.*besselk(1, w_out, 1).*besseli(1, w_in, 1);

    zin = (1i*omega*mu_c/(2*pi)).*(1/w_in).*(N./D);
end

end
