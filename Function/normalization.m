function[new_data]=normalization(old_data)
%��ÿ�������ڲ�ͬ�����е���ֵ�����
[~,n] = size(old_data);
new_data = old_data;

max_element=max(old_data(:,1:n),[],'omitnan');   %�����1�е���N-1�е�ÿ�����ݵ���Сֵ
min_element=min(old_data(:,1:n),[],'omitnan');   %�����1�е���N-1�е�ÿ�����ݵ���Сֵ
for i=1:n
    if max_element(:,i) ~= 0
        new_data(:,i)=(old_data(:,i)-min_element(i))/(max_element(i)-min_element(i));
    end
end

