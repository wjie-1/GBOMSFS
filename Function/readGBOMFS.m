function [temp,selections,time]= readGBOMFS(router, Num ,Smooth, dataset, gaintshold1, gaintshold2, threshold)
warning('off')
addpath(genpath(pwd));


originRoute = router;

data = load([originRoute, dataset '.mat']);  % 假设数据集文件存放在 data 文件夹内


% X = normalize(data.train_data, 'range');
% test_data = normalize(data.test_data, 'range');
% X = normalization(data.train_data); %按列归一化
% test_data = normalization(data.test_data);


X = data.train_data;
test_data = data.test_data;


test_label = data.test_label;
Y = data.train_label;






[n,m] = size(X);
[~,l] = size(Y);



[selections, time] = GBOMFS(X, Y, gaintshold1, gaintshold2, threshold);



fprintf('Selections:\n');
disp(selections); 
disp(time);



Num = 10;Smooth = 1; 
[Prior,PriorN,Cond,CondN]=MLKNN_train(X(:, selections), Y',Num,Smooth);
[HammingLoss,RankingLoss,Coverage,Average_Precision,OneError,macrof1,microf1,Outputs,Pre_Labels]=...
    MLKNN_test(X(:, selections), Y', test_data(:, selections), test_label',Num,Prior,PriorN,Cond,CondN);


% [Average_Precision,Coverage,HammingLoss,OneError,RankingLoss] = mlknn(X(:, selections),Y',test_data(:, selections),test_label');
temp = [Average_Precision, HammingLoss,OneError, RankingLoss, Coverage]; % 五个指标的结果
disp(temp);

%
folderPath = ['.\RES\result_',num2str(threshold),dataset]; % 将 yourFolderName 替换为实际的文件夹名
% 如果文件夹不存在，先创建文件夹
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end

% 拼接完整的保存路径，使用 datasetName 来生成唯一的文件名
savePathData = fullfile(folderPath, [dataset  '_selection_' num2str(gaintshold1) '_' num2str(gaintshold2) '.mat']);
% 使用 save 函数保存结构数组到指定路径
save(savePathData, 'selections','temp','time','threshold');




end



