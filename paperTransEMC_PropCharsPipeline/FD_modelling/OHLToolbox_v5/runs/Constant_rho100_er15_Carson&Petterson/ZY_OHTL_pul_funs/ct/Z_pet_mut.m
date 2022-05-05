% Z_pet_mut.m
function [Zg_pet_mutual]=Z_pet_mut(num,h,d,e_g,m_g,sigma_g,omega)
% This function Calculates the Mutual Earth Impedance by Semlyen's
% Approxiamtion

% input: h1,h2,d,e_g,_m_g,sigma_g,omega
% output: Z_mutual
m0=4*pi*1e-7;
e0=8.854187817e-12;
b=sqrt(1i*omega*m0*sigma_g+omega.^2*(e0*m0-e_g*m0));

          % Mutual Impedance

Zg_pet_mutual=zeros(num,num);

for x=1:1:num
    for y=1:1:num
       if x==y
         Zg_pet_mutual(x,y)=0;
       else
          dp=sqrt((h(1,x)+h(1,y)+2./b).^2+d(x,y).^2);
          d2=sqrt(d(x,y).^2+(h(1,x)+h(1,y)).^2);
          term=log(dp./d2);


          Zg_pet_mutual(x,y)=1i*omega.*m0./(2*pi).*term;
       end
     end
end
