function [V]=calc_Vnode(ord,f,Ybranch,Vo1)
% For details see:
% A. I. Chrysochos, T. A. Papadopoulos, G. K. Papagiannis, "Enhancing the Frequency-Domain Calculation of Transients in Multiconductor Power Transmission Lines," 
% Electric Power Systems Research, vol. 122, 2015, p 56-64, 
% doi.org/10.1016/j.epsr.2014.12.024



V=zeros(2*ord,max(size(f))); % (2*ord x max(size(f)))
Ybranch_dis_f=zeros(2*ord,2*ord); % (2*ord x 2*ord)
for k=1:1:max(size(f))
    
    for o=1:1:2*ord
        Ybranch_dis_f(o,:)=Ybranch(k,(o-1)*2*ord+1:o*2*ord); % Convert Ybranch matric to square matrix (2*ord x 2*ord)
    end
    
    Ired=-Ybranch_dis_f(2:length(Ybranch_dis_f),1)*Vo1(k); % Calculate matrix Ired at the given frequency - ((2*ord)-1 x 1) - No current sources are cnosidered
    Yred=Ybranch_dis_f(2:length(Ybranch_dis_f),2:length(Ybranch_dis_f)); % Calculate Õred matrix at the given frequency - ((2*ord)-1 x (2*ord)-1)
    %Vred=inv(Yred)*Ired; % Calculate Vred matrix at the given frequency -
    %((2*ord)-1 x 1)
    Vred=Yred\Ired;
    Vdis=[Vo1(k);Vred]; % Store node voltages at the given frequency - (2*ord x 1)
    V(:,k)=Vdis; % Store node voltages at the whole frequency spectrum - (2*ord x max(size(f)))
end