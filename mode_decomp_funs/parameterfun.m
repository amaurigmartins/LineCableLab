function [F,J]=parameterfun(x,Q,O,ord)

F=[(Q*x(1:2*ord)-[x(2*ord+1)*O -x(2*ord+2)*O;x(2*ord+2)*O x(2*ord+1)*O]*x(1:2*ord));sum((x(1:ord).^2)-(x(ord+1:2*ord).^2))-1;sum(x(1:ord).*x(ord+1:2*ord))-0];

J_temp=[(S-x(ord+1)*O) -x(1:ord)];