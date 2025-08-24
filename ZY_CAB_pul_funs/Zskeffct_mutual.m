% Zskeffct_mutual.m
function [zm]=Zskeffct_mutual(radius_ex,radius_in,sigma_c,mu_r,omega,varargin)
% Calculation of mutual impedance between the outer and inner surfaces of a
% tubular conductor

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
    cschTerm = csch(m.*(radius_ex-radius_in));
    zm = m./(sigma_c.*pi.*(radius_in+radius_ex)).* cschTerm;
    
else
    w_out=m.*radius_ex;
    w_in=m.*radius_in;
    s_in = exp(abs(real(w_in)) - w_out);
    s_out = exp(abs(real(w_out)) - w_in);
    sc = s_in./s_out; %should be applied to all besseli() involving w_in
    D = besseli(1, w_out, 1).*besselk(1, w_in, 1) - sc.*besselk(1, w_out, 1).*besseli(1, w_in, 1);
    zm = (1./(2.*pi.*radius_ex.*radius_in.*sigma_c.*D.*s_out));
end

end
