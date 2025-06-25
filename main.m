% clear; clc;


warning('off')
addpath(genpath(pwd));


router = 'data\';
k = 10;Smooth = 1;
datasetName={
'Birds';
    }'; 


for num = [1]%:3

    if num == 1
        dataSetNum = length(datasetName);

        folderPath = 'result_GBOMFS2'; % 将 yourFolderName 替换为实际的文件夹名
        % 如果文件夹不存在，先创建文件夹
        if ~exist(folderPath, 'dir')
            mkdir(folderPath);
        end

%         % 创建一个空的结构数组，用于保存数据
%         dataStruct = struct( 'temp1', {}, 'selections', {}, 't', {}, 'gainshold',{});


        %% 特征选择
        for i = 1:dataSetNum

            dataset = datasetName{i};
            Temp=[];
            threshold = 0.8 ;
            for gaintshold1 = 0.001:0.002:0.009
                for gaintshold2 = 0.001:0.002:0.009
                
                [temp,selections,time] = readGBOMFS(router, k, Smooth, dataset, gaintshold1,gaintshold2, threshold);
%                 disp('finished---GBOMFS--')
                Temp = [Temp;temp,gaintshold1,gaintshold2;]
                disp(temp);
                end
            end
                            savePathData = fullfile(folderPath, [dataset  '_selection.mat']);
                % 使用 save 函数保存结构数组到指定路径
%                 save(savePathData, 'dataStruct');
                save(savePathData, 'Temp','time');
            %
            disp(Temp);
                disp('finished---GBOMSFS--')

        end
    end
end