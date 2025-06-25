function [simMatrix,Matrix,ppf] = computeFSimilarity(dataA, dataB, pvalue)
    % 获取数据的行数和列数
    [rows, cols] = size(dataA);
    [nK, ~] = size(dataB);
  

    % 计算每个条件属性的标准差
    sigmas = std(data,"omitnan");

    % 初始化一个空矩阵来存储所有属性的相似性
    simMatrix = zeros(rows,nK)

    % 针对每个条件属性计算相似性矩阵
    for a = 1:cols

        % 创建一个上三角矩阵，不包括对角线
        simMatrixRow = zeros(1, rows * (rows - 1) / 2);
        
        % 计算上三角的相似性 
        idx = 1;
        for x = 1:rows
            values = [dataA(x,:);dataB];
            stdValue = std(values)
            for y = x+1:rows % 不包括x=y的情况
                % 应用提供的公式计算相似性
                value1 = (data(x,a) - data(y,a)) / sigmas(a);
                value2 = (data(y,a) - data(x,a)) / sigmas(a);
                if isnan(data(x,a)) || isnan(data(y,a))
                    simMatrixRow(idx) = 1;
                else
                    simMatrixRow(idx) = max(min(1 - value1, 1 - value2), 0);

                    idx = idx + 1;
                end
            end
        end



        % 自适应阈值设置（均值+0.5标准差）
        mean_sim = mean(simMatrixRow(:));
        std_sim = std(simMatrixRow(:));
        ppf = ( mean_sim + 3 * std_sim)/(mean_sim + 3*0.5);


        % 将相似性行向量放入大矩阵中
        Matrix(a, :) = (simMatrixRow >= ppf);
    end



    % 将相似性矩阵中小于pvalue的值设为0，其余值设为1
    simMatrix = Matrix;
end
