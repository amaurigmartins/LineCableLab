% Z_snd_mut.m
function [Zg_sln_mutual]=Z_sln_mut(num,h,d,e_g,m_g,sigma_g,omega)
% This function Calculates the Mutual Earth Impedance by Semlyen's
% Approxiamtion

% input: h1,h2,d,e_g,_m_g,sigma_g,omega
% output: Z_mutual
m0=4*pi*1e-7;

Zg_sln_mutual=zeros(num,num);

for x=1:num
    for y=1:num
        if x==y
           Zg_sln_mutual(x,y)=0;
        else
            
            
       % Mutual Impedance
       ter1ar_mut=(h(1,x)+h(1,y));
       ter1par_mut=pi*(d(x,y)^2+(h(1,x)+h(1,y))^2);
       ter1_mut=ter1ar_mut./ter1par_mut;

       ter2ar_mut=1i*omega.*m_g;
       ter2par_mut=sigma_g+1i*omega.*e_g;
       ter2_mut=sqrt(ter2ar_mut./ter2par_mut);

       Zg_sln_mutual(x,y)=ter1_mut.*ter2_mut;
       
        end
    end
end
       
