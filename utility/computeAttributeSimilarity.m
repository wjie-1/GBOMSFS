function [colDistances,Dist] = computeAttributeSimilarity(data,S, Dist, radius)

%X是原始数据矩阵，R是已选属性S的矩阵,Dist是距离矩阵
%返回邻域容差相似矩阵，以及距离矩阵
R = data(:,S);
% 获取数据的行数和列数
[n,m] = size(R);  %对象数量,m实际是1

colDistances  = ones(n, n);


for j = 1:n
    for k = j+1:n
        if isnan(R(j,1)) || isnan(R(k,1))
            continue;
        else
            Dist(j,k) = Dist(j,k) + abs(R(j,1) - R(k,1));  % 计算曼哈顿距离
            Dist(k,j) = Dist(j,k);
            %             dist = dist + sqrt(X(j,i) - X(k,i)); %计算欧式距离
            %             dist = dist + max(X(j,i) - X(k,i));
            if (Dist(j,k) - radius)> 1e-5
                colDistances(j,k) = 0;
                colDistances(k,j) = 0;  % 利用对称性
            end
        end
        
    end
    
end

end
