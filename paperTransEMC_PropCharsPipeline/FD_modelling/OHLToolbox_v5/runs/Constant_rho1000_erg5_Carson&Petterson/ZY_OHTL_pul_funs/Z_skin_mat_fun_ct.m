%Internal impedance matrix
function [Zin_mat]=Z_skin_mat_fun_ct(num,rad_ex,rad_in,sigma_w,mrw,omega,Geom)
m0=4*pi*1e-7;        % Henry's/meters

%Create num*num matrix
Zin_mat=zeros(num,num);
Zskin_self=zeros(num,num);
Zins=zeros(num,num);

%% Skin effect
for k=1:1:num
    %Calculate the self Impedance
    Zskin_self(k,k)=skeffct_tb_fun_ct(rad_ex(k),rad_in(k),sigma_w(k),mrw(k),omega);
    if isnan(Zskin_self(k,k))
        Zskin_self(k,k)=0;
    end
    Zin_mat(k,k)=Zskin_self(k,k);
end


for k=1:1:num
    if isnan(Geom(k,8))
        Zins(k,k)=0;
        Zin_mat(k,k)=Zins(k,k)+Zin_mat(k,k);
    else
        Zins(k,k)=1i*omega*m0/(2*pi)*log(Geom(k,8)/Geom(k,5));
        Zin_mat(k,k)=Zins(k,k)+Zin_mat(k,k);
    end   
end


