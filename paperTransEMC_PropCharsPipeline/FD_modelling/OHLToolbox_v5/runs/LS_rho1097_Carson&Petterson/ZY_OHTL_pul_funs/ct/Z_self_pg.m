%Zpg Self function
function [Zpg_self]=Z_self_pg(num,h,d,r,omega)

m0=4*pi*1e-7;
e0=8.854187817e-12;

Zpg_self=zeros(num,num);




for x=1:1:num
    for y=1:1:num
        if x==y
             % Self Impedance
             term_self=log((2*h(1,y))./r(x,y));
             Zpg_self(x,y)=j*omega.*m0./(2*pi).*term_self;
        end
    end
end
