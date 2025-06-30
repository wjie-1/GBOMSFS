function [selections, time] = GBOMFS(X,Y,gaintshold1, gaintshold2, threshold)
tic;

[n,F] = size(X);
[len] = size(Y,2);

YY = (Y==1);

[L, C_combined] = LabelDistribution(YY, 0.5);

%
% 调用 Granular Ball 函数
% [centers, gb_list_final, sample_num] = getGranularBall(L);

[centers, gb_list_final, sample_num] = getGranularBall(L);


GB_num = size(sample_num,2);

% 打印结果
% disp('Centers of Granular Balls:');
% disp(centers);
% 
% disp('\nNumber of Samples in Each Ball:');
% disp(GB_num);
% 
% disp('\nGranular Ball Details:');
% for i = 1:length(gb_list_final)
%     disp(['Ball ', num2str(i), ':']);
%     disp(gb_list_final{i});
% end

[C_matrix, neighbors] = getGB_KNeighbor(centers);

[Q] = getGB_Quantity(gb_list_final, centers, C_matrix, neighbors);
%% online streaming FS



%% 初始化全局变量
S = [];  %已选特征


LastRCI = 0;  %记录总CRCI


fea_RCIList = zeros(1,F); %用于保存所有特征的RCI
GB_Card_Sum = zeros(F,GB_num); %用于保存所有特征的CRCI


%% 在线流特征选择OSFS
% for i = 1:F
    %% 在线流特征选择OSFS
    randkey = 10;
    rng(randkey)
    randFi=randperm(F);
for i = randFi

    gb_card_fi = [];
    for r = 1: GB_num
        [gb_cardi] = getSimilarityGBcard(r, gb_list_final, X(:,i), neighbors, threshold);
        gb_card_fi(end+1) = gb_cardi;
    end

    GB_Card_Sum(i,:) = gb_card_fi;

    feaRCI = FCG(gb_list_final, gb_card_fi, GB_num, n, Q);
    fea_RCIList(i) =  feaRCI;

    %流入特征,计算特征的相似关系矩阵
    if isempty(S)  %若i-1时刻暂无选入的特征
        S(end+1) = i;

        LastRCI = LastRCI + feaRCI;
        MeanRCI = LastRCI / length(S);
        continue;

    else



        %% 强一致增益特征，重要性分析

        if feaRCI > MeanRCI

            S(end+1) = i;

            LastRCI = LastRCI + feaRCI;
            MeanRCI = LastRCI / length(S);


        else


            value11 = [];
            value22 = [];
            for k = S
                value11(end+1) = FCCG(gb_list_final, GB_Card_Sum(i,:), GB_Card_Sum(k,:), GB_num, n, Q);
%                 value22(end+1) = FCCG(gb_list_final, GB_Card_Sum(k,:), GB_Card_Sum(i,:), GB_num, n, Q);
            end
            temp1 = sum(value11)/length(S);


            metric = temp1/feaRCI
            flag = 0;
            [v,idx] = min(value11);
            if metric >= gaintshold1
                [~,imax] = max(value11);
                value2 = FCCG(gb_list_final, GB_Card_Sum(S(imax),:), GB_Card_Sum(i,:), GB_num, n, Q);
                tt = (value11(imax)/feaRCI  - value2/ fea_RCIList(S(imax))) 

                if  tt >= gaintshold2
                    S(idx) = i;

                    LastRCI = LastRCI + feaRCI - fea_RCIList(S(idx));
                    MeanRCI = LastRCI / length(S);


                    flag = 1;
                    continue;
                end

                if flag == 0 
                    S(end+1) = i;

                    LastRCI = LastRCI + feaRCI;
                    MeanRCI = LastRCI / length(S);
                end
            end
         
        end

        %             end

        %             end

    end
end


selections = S;

% end
time = toc;
end



%% calculate K-Neighbor GranularBall
function [C_matrix, neighbors] = getGB_KNeighbor(centers)

C = 1 - pdist(centers,'cosin');
% C = 1 - pdist(centers,'euclidean');

n = size(centers,1);

C_matrix = squareform(C) + eye(n,n); % 转换为方阵

threshold = 0.8;  % 相似度阈值
n = size(C_matrix, 1);  % 粒球数量
neighbors = cell(n, 1);  % 存储每个粒球的近邻粒球

for i = 1:n
    neighbors{i} = find(C_matrix(i, :) >= threshold);  % 找到相似度 >= 0.8 的粒球
    neighbors{i}(neighbors{i} == i) = [];  % 去除自身
end
end

%% calculate GranularBall Quantity
function [Q] = getGB_Quantity(gb_list_final, centers, C_matrix, neighbors)
n = size(centers,1);
Q = [];

for r =1:n
    gbnum = size(gb_list_final{r},1);
    center = centers(r,:);


    % 扩展 x_i 和 gbK 数据
    m1 = repmat(center, gbnum, 1); % 扩展 x_i
    m2 = gb_list_final{r}(:,1:end-1); % gbK 数据

    % 计算欧氏距离
    distances = sum (sqrt(sum((m1 - m2).^2, 2)) );% 对每一行求和


    Q1 = (exp(-distances/gbnum) );
    Q(r) = Q1 ;

end


end

