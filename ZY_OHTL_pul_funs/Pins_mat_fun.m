%Internal impedance matrix
function [Pin_mat]=Pins_mat_fun(num,Geom)
e0=8.854187817e-12;  % Farads/meters
%Create num*num matrix
Pin_mat=zeros(num,num);

for k=1:1:num
    if isnan(Geom(k,8))
        Pin_mat(k,k)=0;
    else
%         Pin_mat(k,k)=log(Geom(k,8)/Geom(k,5))/(2*pi*e0*Geom(k,10));
        Pin_mat(k,k)=log(Geom(k,8)/Geom(k,5))/Geom(k,10);
    end   
end