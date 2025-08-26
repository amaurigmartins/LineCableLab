function [G] = calcGreenFuns(src,tgt,rho)


coords_LT=src;
n_LT=length(src);
coords_DT=tgt;
n_DT=length(tgt);


%basic Green function for homogeneous soil

for i=1:n_LT
    for j=1:n_DT
        xT=coords_LT(i,1);yT=coords_LT(i,2);
        xD=coords_DT(j,1);yD=coords_DT(j,2);
        r=sqrt((xT-xD)^2+(yT-yD)^2);
        G(j,i) = rho/(2*pi*r); %

    end
end


end