% Z_snd_mut.m
function [Zg_der_mutual]=Z_der_mut(num,h,d,e_g,m_g,sigma_g,omega)
% This function Calculates the Mutual Earth Impedance by Semlyen's
% Approxiamtion

% input: num,h,d,e_g,_m_g,sigma_g,omega
% output: Z_mutual
m0=4*pi*1e-7;
p=1./sqrt(j*omega.*m_g*sigma_g);
%p=1./sqrt(j*omega.*m_g.*(sigma_g+j*omega.*e_g));

Zg_der_mutual=zeros(num,num);
% Mutual Impedance

for x=1:num
    for y=x+1:num
        if x==y
            Zg_der_mutual(x,y)=0;
        else
            terar_mut=sqrt((h(1,x)+h(1,y)+2*p).^2+d(x,y)^2);

            terpar_mut=sqrt((h(1,x)+h(1,y))^2+d(x,y)^2);

            Zg_der_mutual(x,y)=j*omega.*m_g./(2*pi).*log(terar_mut./terpar_mut);
            Zg_der_mutual(y,x)=Zg_der_mutual(x,y);
         
        end
     end
end
