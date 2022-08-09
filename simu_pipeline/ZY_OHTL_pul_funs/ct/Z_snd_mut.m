% Z_snd_mut.m
function [Zg_mutual]=Z_snd_mut(num,h,d,e_g,m_g,sigma_g,omega)
% This function Calculates the Mutual Earth Impedance by Sunde's
% Approxiamtion

% input: h1,h2,d,e_g,_m_g,sigma_g,omega
% output: Z_mutual
m0=4*pi*1e-7;
gamma_g=sqrt(1i*omega.*m_g.*(sigma_g+1i*omega.*e_g));

% Mutual Impedance
Zg_mutual=zeros(num,num);


for x=1:1:num
    for y=1:1:num
        if x==y 
            
            ar_mut=1;
            par_mut=1;
            %Zg mutual for x=y is the self impedance
            %We set it 0 in the mutual matrix
            %log1=0
            Zg_mutual(x,y)=1i*omega.*m0/(4*pi).*log(ar_mut./par_mut);
        else
            

            ar_mut=(1+gamma_g.*((h(1,x)+h(1,y))/2)).^2+(gamma_g.*(d(x,y)/2)).^2;
            par_mut=(gamma_g.*((h(1,x)+h(1,y))/2)).^2+(gamma_g.*(d(x,y)/2)).^2;

            Zg_mutual(x,y)=1i*omega.*m0/(4*pi).*log(ar_mut./par_mut);
            Zg_mutual(y,x)=Zg_mutual(x,y);
            

        end
    end

end

