function Y=add_branch(n1,n2,y,Y);
%
% Y=add_branch(n1,n2,y,Y);
%
% The function adds the contribution fom an admittance y between nodes 
% n1 and n2 to the admittance matrix Y. Node number 0 denotes ground.
%
%  Download site:
%  http://www.energy.sintef.no/Produkt/VECTFIT/index.asp 
%
%  17.03.2013. Bjorn Gustavsen, SINTEF Energy Research, Norway.

if n1==0 & n2==0
  return
end
if n1==0
  Y(n2,n2)=Y(n2,n2)+y;
elseif n2==0  
  Y(n1,n1)=Y(n1,n1)+y;
else  
  Y(n1,n1)=Y(n1,n1)+y;
  Y(n2,n2)=Y(n2,n2)+y;
  Y(n1,n2)=Y(n1,n2)-y;
  Y(n2,n1)=Y(n2,n1)-y;
end  