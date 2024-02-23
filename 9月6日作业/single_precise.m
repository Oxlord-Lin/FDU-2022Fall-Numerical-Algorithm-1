k=10.^8;
Sum1=0;
Sum2=0;
for n=1:k;
    a=single(1/(n));
    Sum1=Sum1 + a;
    if Sum1==Sum2;
        End_of_n=n-1;%记录最后一个可以表示的a的下标n-1
        break %如果Sum1==Sum2，即1/n在单精度计算下，发生了下溢，变成了0，就停止循环
    end
    Sum2=Sum1;
end
Sum1
n-1