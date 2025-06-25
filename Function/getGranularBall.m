

%% 粒球划分函数
function [centers, gb_list_final, sample_num] = getGranularBall(data)
    [n, m] = size(data);
    
    % 添加索引列
    index = (1:n)';  % 创建索引列，MATLAB中索引从1开始
    data_index = [data, index];  % 将数据和索引列拼接
    
    gb_list_temp = {data_index};  % 使用cell数组存储数据
    if sqrt(n) > 64
        min_split_num = 64;
    elseif sqrt(n) > 8
        min_split_num = 8;
    else
        min_split_num = round(sqrt(n));
    end

    % 数据分割（递归分割直到不再变化）
    while true
        ball_number_old = length(gb_list_temp);
        gb_list_temp = division_ball(gb_list_temp, min_split_num);
        ball_number_new = length(gb_list_temp);
        if ball_number_new == ball_number_old
            break;
        end
    end

    % 计算每个球的半径
    radius = [];
    for i = 1:length(gb_list_temp)
        if size(gb_list_temp{i}, 1) >= 2
            radius = [radius, get_radius(gb_list_temp{i})];  % 存储每个球的半径
        end
    end

    % 计算半径的中位数和均值
    radius_median = median(radius);
    radius_mean = mean(radius);
    radius_detect = max(radius_median, radius_mean);

    % 标准化球的大小
    while true
        ball_number_old = length(gb_list_temp);
        gb_list_temp = normalized_ball(gb_list_temp, radius_detect);
        ball_number_new = length(gb_list_temp);
        if ball_number_new == ball_number_old
            break;
        end
    end

    gb_list_final = gb_list_temp;

    % 获取每个球的样本数量和球心
    sample_num = cellfun(@(x) size(x, 1), gb_list_final);
    gb_centers = cellfun(@(x) mean(x(:, 1:end-1), 1), gb_list_final, 'UniformOutput', false);
%     centers = cell2mat(gb_centers)';
    centers = vertcat(gb_centers{:});  % 按行拼接
end

% 计算球的半径
function radius = get_radius(gb)
    [n, m] = size(gb);
    gb = gb(:, 1:m-1);  % 去除最后一列
    center = mean(gb, 1);  % 计算球的中心
    diffMat = bsxfun(@minus, gb, center);  % 计算每个点到中心的差值
    sqDistances = sum(diffMat.^2, 2);  % 计算平方距离
    distances = sqrt(sqDistances);  % 计算欧几里得距离
    radius = max(distances);  % 半径为最大距离
end

% 分割球：使用k-means方法
function [ball1, ball2] = spilt_ball(data)
    [n, m] = size(data);
    cluster = kmeans(data(:, 1:m-1), 2, 'Start', 'plus', 'Replicates', 5);
    ball1 = data(cluster == 1, :);
    ball2 = data(cluster == 2, :);
%     gb_list_new = {ball1, ball2};
end

% 分割球：基于最大距离
function [ball1, ball2] = spilt_ball_2(data)
    [n, m] = size(data);
    X = data(:, 1:m-1);  % 去除最后一列
    D = pdist2(X, X, 'euclidean');  % 计算距离矩阵
    [r, c] = find(D == max(D(:)));  % 找到最大距离的索引
    r1 = r(2);  % 第二个最大值的行
    c1 = c(2);  % 第二个最大值的列
    
    cluster = zeros(n, 1);
    for j = 1:n
        if D(j, r1) < D(j, c1)
            cluster(j) = 1;
        else
            cluster(j) = 2;
        end
    end
    ball1 = data(cluster == 1, :);
    ball2 = data(cluster == 2, :);
%     gb_list_new = {ball1, ball2};
end

% 计算每个球的稠密度体积
function result = get_density_volume(gb)
    gb = gb(:, 1:end-1);  % 去除最后一列
    num = size(gb, 1);
    center = mean(gb, 1);  % 计算中心
    diffMat = bsxfun(@minus, gb, center);  % 计算每个点与中心的差值
    sqDistances = sum(diffMat.^2, 2);  % 计算平方距离
    distances = sqrt(sqDistances);  % 计算欧几里得距离
    sum_radius = sum(distances);  % 总半径

    result = (1 + exp(-sum_radius/num) )^-1;
%     result = sum_radius;  % 结果为总半径，如果为0，返回样本数
end

% 分割球：根据稠密度进行分割
function gb_list_new = division_ball(gb_list, n)
    gb_list_new = {};
    for i = 1:length(gb_list)
        gb = gb_list{i};
        if size(gb, 1) >= n
            [ball_1, ball_2] = spilt_ball_2(gb);
            density_parent = get_density_volume(gb);
            density_child_1 = get_density_volume(ball_1);
            density_child_2 = get_density_volume(ball_2);
            w = size(ball_1, 1) + size(ball_2, 1);
            w1 = size(ball_1, 1) / w;
            w2 = size(ball_2, 1) / w;
            w_child = (w1 * density_child_1 + w2 * density_child_2);  % 加权稠密度
            
            if w_child < density_parent
                gb_list_new = [gb_list_new, {ball_1, ball_2}];
            else
                gb_list_new = [gb_list_new, {gb}];
            end
        else
            gb_list_new = [gb_list_new, {gb}];
        end
    end
end

% 标准化球
function gb_list_temp = normalized_ball(gb_list, radius_detect)
    gb_list_temp = {};
    for i = 1:length(gb_list)
        gb = gb_list{i};
        if size(gb, 1) < 2
            gb_list_temp = [gb_list_temp, {gb}];
        else
            [ball_1, ball_2] = spilt_ball_2(gb);
            if get_radius(gb) <= 2 * radius_detect
                gb_list_temp = [gb_list_temp, {gb}];
            else
                gb_list_temp = [gb_list_temp, {ball_1, ball_2}];
            end
        end
    end
end
