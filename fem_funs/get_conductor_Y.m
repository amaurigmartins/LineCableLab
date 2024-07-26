function [out] = get_conductor_Y(Geom,j,id)

if nargin == 2
    compute_integral=true;
else
    compute_integral=false;
end

global HandleToFEMM;

x_c = Geom(j,2);
y_c = Geom(j,3);
r_ext = Geom(j,5);
r_ins = Geom(j,8);

DELTA=1e-3;
r_contour=max([r_ext r_ins])+DELTA;
co_seteditmode('contour');
co_clearcontour();
window_size=3*r_contour;
co_zoom(x_c-window_size,y_c-window_size,x_c+window_size,y_c+window_size);

if compute_integral
    theta=[0:(pi/200):2*pi];
    
    for jj=1:length(theta)
        x=x_c+r_contour*cos(theta(jj));
        y=y_c+r_contour*sin(theta(jj));
        co_addcontour(x,y);
        % pause(0.01);
    end
    
    vals=co_lineintegral(1);
    out = -vals(1);
    
else
    vals = co_getconductorproperties(id);
    out = vals(2);
end

end
