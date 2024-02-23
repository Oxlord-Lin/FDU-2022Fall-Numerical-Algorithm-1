c=1;
[a,b] = block_size(temp)

function [l,m] = block_size(T)
n = size(T);
n = n(1,1);
% 分块【分块部分有bug】
    % 首先确定m
    flg = 0;
    for i = n-1:-1:1
        T(i+1,i) 
        if T(i+1,i) ~= 0
            flg = i; %flg所指的是T22的最右侧元素
            break
        end
    end
    if flg ~=0
        m = n-flg-1;
        %only for test
    else
        m = n;
        %only for test
        % 注意：此时已经变成对角阵，应当结束循环，跳出循环输出有关信息
    end
    
    % 确定l【这一步有bug？】
    flg = 0;
    for j = n-m-1:-1:1 
        % 这一步是在确定m后且没有结束循环的情况下进行的，也就是仍有一个三对角阵T22
        % 此时n-m所指向的就是三对角阵最右侧的位置
        if T(j+1,j) == 0
            flg  = j; % 此时j所指的是T11最右侧的位置
            break;
        end
    end
    if flg == 0
        % 此时从1到i+1都是三对角阵的范围
        l = 0;
    else
        % 此时从1到j都是T11的部分
        l = flg;
    end
end