function [V]=calc_Vnode_sheath_inj(ord,f,Ybranch,Vo1)

V=zeros(2*ord,max(size(f))); % (2*ord x max(size(f)))
Ybranch_dis_f=zeros(2*ord,2*ord); % (2*ord x 2*ord)

for k=1:1:max(size(f))
    
    for o=1:1:2*ord
        Ybranch_dis_f(o,:)=Ybranch(k,(o-1)*2*ord+1:o*2*ord);
    end
    Y_temp=Ybranch_dis_f(:,2);
    Y_temp(2)=[];
        
    Ired=-Y_temp*Vo1(k);
    
    Yred=Ybranch_dis_f;
    Yred(2,:)=[];
    Yred(:,2)=[];
    
    Vred=Yred\Ired;
    Vdis=[Vred(1);Vo1(k);Vred(2:(2*ord)-1)];
    V(:,k)=Vdis;
end