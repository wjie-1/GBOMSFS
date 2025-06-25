function [label_distributions, C_combined] = LabelDistribution(Y, alpha)
    % 输入:
    %   Y: n x m 的多标签矩阵 (n 样本数, m 标签数)
    %   alpha: 权重参数 (0 <= alpha <= 1)，用于平衡余弦相似度和共现信息
    % 输出:
    %   label_distributions: n x m 的标签分布矩阵
    %   C_combined: m x m 的综合标签相关性矩阵

    [n_samples, n_labels] = size(Y);

    % 计算余弦相似度
    cosine_sim = 1 - pdist(Y', 'cosine');
    cosine_sim = squareform(cosine_sim) + eye(n_labels, n_labels);

    % 计算共现矩阵
    C_cooccurrence = Y' * Y;
    C_cooccurrence = C_cooccurrence - diag(diag(C_cooccurrence)); % 去除对角线元素
    max(C_cooccurrence(:))
    C_cooccurrence_norm = C_cooccurrence / max(C_cooccurrence(:)); % 归一化共现矩阵

    % 结合余弦相似度和共现信息
    C_combined = alpha * cosine_sim + (1 - alpha) * C_cooccurrence_norm;

    % 处理全零标签
    zero_labels = all(Y == 0, 1); % 找到全零标签
    C_combined(zero_labels, :) = 0; % 将这些标签的相似度设为 0
    C_combined(:, zero_labels) = 0; % 对称处理

    
    P_positive = sum(Y,1)/n_samples;
    P_negative = (n_samples - sum(Y,1))/n_samples;
    
    % 为每个样本生成标签分布
    label_distributions = zeros(n_samples, n_labels);
    for i = 1:n_samples
        weights = zeros(1, n_labels);
        for j = 1:n_labels
            if Y(i, j) == 1 % 如果样本有标签 j
                weights = weights + C_combined(j, :); % 增加与标签 j 相关的标签权重
            end
        end
%         WP = weights;
%         for j = 1:n_labels
%             if Y(i, j) == 1 % 如果样本有标签 j
%                 WP(j) = weights(j)*P_positive(j); % 增加与标签 j 相关的标签权重
%             else
%                 WP(j) = weights(j)*P_negative(j);
%             end
%         end
%    
%         if sum(WP) > 0 % 避免除以零
%             WP = WP / sum(WP); % 归一化
%         end
%         label_distributions(i, :) = WP;

        if sum(weights) > 0 % 避免除以零
            weights = weights / sum(weights); % 归一化
        end
        label_distributions(i, :) = weights;
    end

    % 输出标签相关性矩阵和生成的标签分布
%     disp('综合标签相关性矩阵:');
%     disp(C_combined);
%     disp('生成的标签分布:');
%     disp(label_distributions);
end