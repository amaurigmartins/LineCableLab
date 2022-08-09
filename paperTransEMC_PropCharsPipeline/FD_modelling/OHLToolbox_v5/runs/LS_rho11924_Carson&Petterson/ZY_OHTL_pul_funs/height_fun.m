% Set the conductors at the same height 

function [h]=height_fun(num,Geom)
h=zeros(1,num);

for k=1:1:num
    h(1,k)=Geom(k,3);
end

