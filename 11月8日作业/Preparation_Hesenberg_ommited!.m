%可以用于获得矩阵A的上海森伯格形式A_hes，并且获得Householder反射子House
%二者满足的关系是：A_hess=House*A*House
%Householder是对称正交阵
[B,H]=hessen(A)
B
H
H*A*H
function [A,House_all]=hessen(A)
[m,~]=size(A);
House_all=eye(m,m);
for k=1:m-2
    v=A(k+1:m,k);
    v=v+sign(v(1,1))*norm(v);
    v=v/norm(v);
    house_local=eye(m-k)-2*(v*v');
    house=[eye(k),zeros(k,m-k);zeros(m-k,k),house_local];
    A=house*A*house;
    House_all=house*House_all;
end
end





    