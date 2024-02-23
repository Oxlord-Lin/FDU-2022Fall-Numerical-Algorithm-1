%% 生成A
m = 500;
n = 300;
Q=zeros(m,n);
R=zeros(n,n);
A=randn(m,n)+(sqrt(-1))*randn(m,n);
A0 = A;
%% 用G-S方法进行QR分解
for j = 1:n
    for i = 1:j-1
        R(i,j)=Q(1:m,i)'*A(1:m,j);
        A(1:m,j)=A(1:m,j)-Q(1:m,i).*R(i,j); %改进的G-S方法
    end
    R(j,j)=norm(A(1:m,j));
    Q(1:m,j)=A(1:m,j)./R(j,j);
end
% imagesc(abs(Q'*Q-eye(n)));
% colorbar
norm(abs(Q'*Q-eye(n)))
%% 使用cholesky方法进行QR分解
% R=chol(A'*A);
% Q=A*inv(R);
% %imagesc(abs(Q'*Q-eye(n)));
% colorbar
% 自己写Cholesky
A = A0;
S = A'*A;
S0 = S; % 保留原始数据
for i = 1:n
    for j = i+1:n
        S(j,j:n) = S(j,j:n) - S(i,j:n)*S(i,j)'/S(i,i); % 此处勿忘转置S（i，j）
    end
    S(i,i:n) = S(i,i:n)./sqrt(S(i,i));
end
for i = 1:n
    for j = 1:i-1
        S(i,j) = 0;
    end
end
R = S;
Q = A/R;
imagesc(abs(Q'*Q-eye(n)));
colorbar;
norm(abs(Q'*Q-eye(n)))
%% 使用householder方法进行QR分解
A2 = A0;
for j =1:n
    v=A2(j:m,j);
    cos=real(v(1,1))/abs(v(1,1));%计算相位
    sin=imag(v(1,1))/abs(v(1,1));%计算相位
    len=sqrt(v'*v);
    if real(v(1,1))>=0%此处相当于对v与一个长度为||v||且只有第一个分量不为零的向量作差
        v(1,1)=v(1,1)+len;%避免舍入误差
    else
        v(1,1)=v(1,1)-len;%避免舍入误差
    end
    H=eye(m-j+1)-2*(v*v')./(v'*v);%基本上找到处理第j列的Householder反射子，但还要调整相位，以保证相乘为实数
    H=H*(cos-sqrt(-1)*sin);%进行相位的修正
    A2(j:m,j:n)=H*A2(j:m,j:n);%Householder反射子只需要作用于schur补即可
end

R=A2;
Q=A/R;
% imagesc(abs(Q'*Q-eye(m)));
% colorbar

% 疑问：使用Householder的方法反而误差非常巨大。
% Householder方法得到的|Q'*Q-eye(n)|大概是G-S和cholesky的10^14倍左右？为什么呢？

%% 复习时重写了一个Householder Triangularization
A = A0;
vector_store = zeros(m,n); % 用于存储各个步骤的反射向量vk
Q = eye(m);%用于储存Householder阵
for k =1:n
    x = A(k:m,k);
    phase_conj = x(1,1)'./abs(x(1,1)); % 将用于调整相位，使得R的对角元为实数
    x(1,1) = x(1,1) + sign(x(1,1))*norm(x);
    vk = x./norm(x);
    vector_store(k:m,k) = vk;
    A(k:m,k:n) = phase_conj*(A(k:m,k:n) - 2*vk*(vk'*A(k:m,k:n)));
    Q(k:m,:) = phase_conj*(Q(k:m,:) - 2*vk*(vk'*Q(k:m,:))); %注意此处要对Q的k到m行都进行行变换
end
R = A(1:n,1:n); % reduced QR factorization

% 如果A是实数阵，那么可以用下述方法计算Q，更加节省空间
% Q = eye(m);
% Q = Q(1:m,1:n);
% for k = n:-1:1 % 计算Q的各列
%     vk = vector_store(k:m,k);
%     Q(k:m,:) = Q(k:m,:) - 2*vk*(vk'*Q(k:m,:));
% end
Q = Q';
Q = Q(:,1:n); % reduced QR factorization
imagesc(abs(Q'*Q-eye(n)));
colorbar;
norm(abs(Q'*Q-eye(n)))
