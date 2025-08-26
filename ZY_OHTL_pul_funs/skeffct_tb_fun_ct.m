function [zin]=skeffct_tb_fun_ct(radius_ex,radius_in,sigma_w,mrw,omega)
%skeffct_tb_fun_ct Calculates the internal impedance of a tubular conductor
%   due to the skin effect.
%
%   Syntax:
%       [zin] = skeffct_tb_fun_ct(radius_ex, radius_in, sigma_w, mrw, omega)
%
%   Inputs:
%       radius_ex - External radius of the conductor (r2) [m]
%       radius_in - Internal radius of the conductor (r1) [m]. If 0, a
%                   solid conductor is assumed.
%       sigma_w   - Electrical conductivity of the conductor (sigma_c1) [S/m]
%       mrw       - Relative magnetic permeability of the conductor (mu_c1)
%       omega     - Angular frequency (omega) [rad/s]
%
%   Output:
%       zin - Internal impedance per unit length [Ohm/m]

% --- Constants
mu0 = 4*pi*1e-7; % Permeability of free space (H/m)
j = 1i;          % Imaginary unit

% --- Implementation

% Common term from equations 2.20 and 2.21
m = sqrt(j * omega * mu0 * mrw * sigma_w);

% w2 is used in both equations
w2 = radius_ex * m;

if radius_in > 0
    % This is a tubular conductor

    % Calculate w1 and s_c1 for the tubular case
    w1 = radius_in * m;
    sc1 = exp(abs(real(w1)) - abs(real(w2)) + w1 - w2);

    % Numerator
    numerator = besseli(0, w2) * besselk(1, w1) + sc1 * besseli(1, w1) * besselk(0, w2);

    % Denominator
    denominator = besseli(1, w2) * besselk(1, w1) - sc1 * besseli(1, w1) * besselk(1, w2);

    % Full expression for zin
    zin = (j * omega * mu0 * mrw) / (2 * pi * w2) * (numerator / denominator);

else
    % This is a solid conductor (radius_in = 0)

    % Numerator
    numerator = besseli(0, w2);

    % Denominator
    denominator = besseli(1, w2);

    % Full expression for zin
    zin = (j * omega * mu0 * mrw) / (2 * pi * w2) * (numerator / denominator);

end

end
