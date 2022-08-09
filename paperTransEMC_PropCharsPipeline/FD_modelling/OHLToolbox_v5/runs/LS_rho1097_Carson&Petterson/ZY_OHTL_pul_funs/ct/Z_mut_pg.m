%Zpg Mutual function
function [Zpg_mut]=Z_mut_pg(num,h,d,r,omega)

m0=4*pi*1e-7;
e0=8.854187817e-12;

Zpg_mut=zeros(num,num);




for x=1:1:num
    for y=1:1:num
        if x~=y
             % Mutual Impedance
             d1=sqrt(d(x,y).^2+((h(1,x)-h(1,y)).^2));
             d2=sqrt(d(x,y).^2+((h(1,x)+h(1,y)).^2));
             term_mut=log(d2./d1);
             Zpg_mut(x,y)=j*omega.*m0./(2*pi).*term_mut;
        else
            Zpg_mut(x,y)=0;
        end
    end
end