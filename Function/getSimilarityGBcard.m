function [gb_card] = getSimilarityGBcard(r, gb_list_final, data, neighbors, threshold)
gbnum = size(gb_list_final{r}, 1);
gbX = gb_list_final{r}(:, end)';
gbX_data = data(gbX, :); % gbX 的样本数据

if size(gbX,1) > 1
value1 = getGBcard(gbX_data, threshold);
else
    value1 = 1;
end
gbK_data = [];
flag = 0;
for j = [neighbors{r}]
    gbK = gb_list_final{j}(:, end)';
    gbK_data = [gbK_data; data(gbK, :)];
    flag = 1;                                             
end

Kcard = 0;
if flag == 1
    stdValue = std([gbX_data; gbK_data]); % 计算标准差
    for i = 1:gbnum
        x_i = gbX_data(i, :); % 当前样本 x_i
        m1 = repmat(x_i, size(gbK_data, 1), 1); % 扩展 x_i
        m2 = gbK_data; % gbK 数据

        % 初始化多属性相似度矩阵
        num_attributes = size(gbK_data, 2); % 属性数量
        C_matrix_final = ones(size(m1, 1), 1); % 初始化为全 1 矩阵
        % 逐属性计算相似度矩阵
        for attr = 1:num_attributes
            % 提取当前属性的数据
            m1_attr = m1(:, attr); % 当前属性的 x_i 数据
            m2_attr = m2(:, attr); % 当前属性的 gbK 数据

            % 计算相似度矩阵
            %             C_matrix = max(min((m2 - (m1 - stdValue)) / stdValue, ((m1 + stdValue) - m2) / stdValue), 0);

            % 计算当前属性的相似度矩阵
            C_matrix_attr = max(min((m2_attr - (m1_attr - stdValue(attr))) / stdValue(attr), ...
                ((m1_attr + stdValue(attr)) - m2_attr) / stdValue(attr)), 0);

            % 对每个属性的相似度矩阵取最小值
            C_matrix_final = min(C_matrix_final, C_matrix_attr);
        end
        C_matrix = C_matrix_final;


        above_threshold = C_matrix(C_matrix >= threshold);                                     
        sum_above = sum(above_threshold); % 大于阈值的值直接相加
        count_above = numel(above_threshold); % 大于阈值的值的个数
%         avg_above = sum_above / max(log2(count_above), 1); % 大于阈值的值求平均，避免除以零
        if sum_above == 0
            avg_above = 0;
        else
            avg_above = sum_above /( threshold*(count_above));
        end
        below_threshold = C_matrix(C_matrix < threshold);
        sum_below = sum(below_threshold); % 小于阈值的值求和
        count_below = numel(below_threshold); % 小于阈值的值的个数
%         avg_below = sum_below / max(log2(count_below), 1); % 小于阈值的值求平均，避免除以零
        
        if sum_below == 0
            avg_below = 0;
        else
            avg_below = sum_below / (threshold*(count_below));
        end

        similarity_values = sum_above + avg_below;
%         similarity_values = avg_above% + avg_below;
        Kcard = Kcard + similarity_values;           
    end
%     Kcard = Kcard/length(neighbors{r});
end                                     
gb_card = value1 + Kcard;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
end

function [gb_card] = getGBcard(gb, threshold)
[C_matrix] = calculateSimilarity(gb, 3);

above_threshold = C_matrix(C_matrix >= threshold);
sum_above = sum(above_threshold); % 大于阈值的值直接相加
count_above = numel(above_threshold); % 大于阈值的值的个数
avg_above = sum_above / max(log2(count_above), 1); % 大于阈值的值求平均，避免除以零

below_threshold = C_matrix(C_matrix < threshold);
sum_below = sum(below_threshold); % 小于阈值的值求和
count_below = numel(below_threshold); % 小于阈值的值的个数
% avg_below = sum_below / max(log2(count_below), 1); % 小于阈值的值求平均，避免除以零
% avg_below = sum_below / max(log2(count_below), 1); 
if sum_below == 0
    avg_below = 0;
else
    avg_below = sum_below / (threshold*(count_below));
end

% avg_below = sum_below; % 小于阈值的值求平均，避免除以零
% gb_card = count_above + sum_below;
gb_card = count_above + avg_below;
end