function [Average_Precision,Coverage,HammingLoss,OneError,RankingLoss]=mlknn(train_data,train_target,test_data,test_target)
Num=10;
Smooth=1;
[Prior,PriorN,Cond,CondN]=MLKNN_train1(train_data,train_target,Num,Smooth);
[HammingLoss,RankingLoss,OneError,Coverage,Average_Precision,Outputs,Pre_Labels]=MLKNN_test1(train_data,train_target,test_data,test_target,Num,Prior,PriorN,Cond,CondN);
