% Z_self_mut_pg.m
function [Zpg_mat]=Z_self_mut_pg(h,d,r,omega,con)
% Series impedances of perfect ground


m0=4*pi*1e-7;
e0=8.854187817e-12;


Zpg_mat=zeros(con,con);


for x=1:con
    for y=1:con
        if x==y
       
            % Self Impedance
            term_self=log((2*h(1,y))./r(x));
            Zpg_mat(x,y)=1i*omega.*m0./(2*pi).*term_self;
       
        else

            % Mutual Impedance
            d1=sqrt(d(x,y).^2+(h(1,x)-h(1,y)).^2);
            d2=sqrt(d(x,y).^2+(h(1,x)+h(1,y)).^2);
            term_mut=log(d2./d1);
            Zpg_mat(x,y)=1i*omega.*m0./(2*pi).*term_mut;

        end
    end
end
