% Z_pet_mut.m
function [Zpg_mat]=Z_self_mut_pg(num,h,d,r,omega)
% This function Calculates the Mutual Earth Impedance by Semlyen's
% Approxiamtion

% input: h1,h2,h3,d,e_g,_m_g,sigma_g,omega
% output: Z_mutual
m0=4*pi*1e-7;
e0=8.854187817e-12;

Zpg_mat=zeros(num,num);




for x=1:num
    for y=1:num
        if x==y
             % Self Impedance
             term_self=log((2*h(1,y))./r);
             Zpg_mat(x,y)=j*omega.*m0./(2*pi).*term_self;

        else
              % Mutual Impedance
             d1=sqrt(d(x,y).^2+((h(1,y)-h(1,y)).^2));
             d2=sqrt(d(x,y).^2+((h(1,y)+h(1,y)).^2));
             term_mut=log(d2./d1);
             Zpg_mat(x,y)=j*omega.*m0./(2*pi).*term_mut;
            
       

        end
    end
end

