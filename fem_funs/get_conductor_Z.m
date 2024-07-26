function [out] = get_conductor_Z(Geom,j,w, id)

if nargin == 3
    compute_integral=true;
else
    compute_integral=false;
end

global HandleToFEMM;

    x_c = Geom(j,2);
    y_c = Geom(j,3);
    r_in = Geom(j,4);
    r_ext = Geom(j,5);
    sig_c=1/Geom(j,6)/1e6;
    r_ins = Geom(j,8);
    r_total=max([r_ext r_ins]);
    dx=(r_ext+r_in)/2;
    window_size=3*r_total;
    mo_zoom(x_c-window_size,y_c-window_size,x_c+window_size,y_c+window_size);
    
if compute_integral
    mo_clearblock();
    mo_selectblock(x_c-dx,y_c);
    S = mo_blockintegral(5);
    A = mo_blockintegral(1)/S; % this is the INTEGRAL across the surface S, ie the current
    J = 1i*w*sig_c*A;
    out = J/sig_c; % yeah, redundant but we gotta respect physics
    
    % Type    Definition
    % 0       A · J
    % 1       A
    % 2       Magnetic field energy
    % 3       Hysteresis and/or lamination losses
    % 4       Resistive losses
    % 5       Block cross-section area
    % 6       Total losses
    % 7       Total current
    % 8       Integral of Bx (or Br) over block
    % 9       Integral of By (or Bz) over block
    % 10      Block volume
    % 11      x (or r) part of steady-state Lorentz force
    % 12      y (or z) part of steady-state Lorentz force
    % 13      x (or r) part of 2× Lorentz force
    % 14      y (or z) part of 2× Lorentz force
    % 15      Steady-state Lorentz torque
    % 16      2× component of Lorentz torque
    % 17      Magnetic field coenergy
    % 18      x (or r) part of steady-state weighted stress tensor force
    % 19      y (or z) part of steady-state weighted stress tensor force
    % 20      x (or r) part of 2× weighted stress tensor force
    % 21      y (or z) part of 2× weighted stress tensor force
    % 22      Steady-state weighted stress tensor torque
    % 23      2× component of weighted stress tensor torque
    % 24      R2 (i.e. moment of inertia / density)
    % 25      x (or r) part of 1× weighted stress tensor force
    % 26      y (or z) part of 1× weighted stress tensor force
    % 27      1× component of weighted stress tensor torque
    % 28      x (or r) part of 1× Lorentz force
    % 29      y (or z) part of 1× Lorentz force
    % 30      1× component of Lorentz torque
else
    vals = mo_getcircuitproperties(id);
    out=vals(2); %2 is the voltage drop which is equal to J/sigma
end
