function [] = inspect_femm_cond(fname,Geom,j)

global HandleToFEMM;
openfemm;
opendocument(fname);

hand1=HandleToFEMM;

[~,~,fext]=fileparts(fname);

if strcmp(fext,'.fec')
    prefix='c';
elseif strcmp(fext,'.fem')
    prefix='m';
else
    error('Don Corleone, I don''t know what to do, I don''t know what to doooo! :(')
end

eval(sprintf('%si_loadsolution();',prefix));
eval(sprintf('%so_showmesh();',prefix));
eval(sprintf('%so_shownames();',prefix));

%Inspect conductors
    x_c = Geom(j,2);
    y_c = Geom(j,3);
    r_ext = Geom(j,5);
    r_ins = Geom(j,8);
    DELTA=1e-5;
    r_contour=sum([r_ext r_ins],'omitnan')+DELTA;
    eval(sprintf('%so_seteditmode(''contour'');',prefix));
    eval(sprintf('%so_clearcontour();',prefix));
    window_size=1*r_contour;
    co_zoom(x_c-window_size,y_c-window_size,x_c+window_size,y_c+window_size);
    theta=[0:(pi/200):2*pi];
    for jj=1:length(theta)
        x=x_c+r_contour*cos(theta(jj));
        y=y_c+r_contour*sin(theta(jj));
        eval(sprintf('%so_addcontour(x,y);',prefix));
        %pause(0.01);
    end

end

