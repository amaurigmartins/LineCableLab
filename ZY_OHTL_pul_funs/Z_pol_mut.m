function [Zg_mutual]=Z_pol_mut(h,d,sigma1,f,con)

% Function for the Mutual Earth Impedance with overhead-underground
% arrangement using Pollaczek formula

% Inputs
% h       : depth of conductor [m]
% d       : distances between conductors [m]
% sigma1  : conductivity of earth [S/m]
% f       : frequency [Hz]

% Output
% out: Mutual earth impedance between overhead/underground conductors [Ohm/m]

% Constants
mu0=4*pi*1e-7;
w=2*pi*f;
m=sqrt(1i.*w.*mu0.*sigma1);

Zg_mutual=zeros(con,con);

% Mutual Impedance
for x=1:con
    for y=x+1:con
        if x~=y
            h1=max([h(1,x) h(1,y)]);
            h2=min([h(1,x) h(1,y)]);
            if sign(h1)~=sign(h2)
                yy=d(x,y);
                F=@(s) exp(-h1.*abs(s)+h2.*sqrt(s.^2+m.^2))./(sqrt(s.^2+m.^2)+abs(s));%to keep consistent with Ametani's paper  10.1109/TEMC.2009.2019953
                zfun=@(s) sum([0 F(s)*cos(s*yy)],'omitnan');
                Jm=integral(zfun,0,Inf,'ArrayValued',true);
                Zg_mutual(x,y)=(1i*w*mu0/pi)*Jm;
                Zg_mutual(y,x)=Zg_mutual(x,y);
            end

        end
    end
end
end

