%Intermediate distance between consecutive conductors

function [d]=distance_fun(num,Geom)

%create a matrix where we will save the distance between the different
%conductors

d=zeros(num,num);

int_cond=1;

for x=1:1:num

    for y=1:1:num
        if x==y
            d(x,y)=0;
            int_cond=1;

        elseif abs(x-y)==int_cond
            d(x,y)=abs(Geom(x,2)-Geom(y,2));
            d(y,x)=d(x,y);
            int_cond=int_cond+1;
        

        end

    end
            int_cond=1;

end