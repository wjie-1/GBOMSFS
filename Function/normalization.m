function[new_data]=normalization(old_data)
%将每个属性在不同论域中的最值差求出
[~,n] = size(old_data);
new_data = old_data;

max_element=max(old_data(:,1:n),[],'omitnan');   %求出第1列到第N-1列的每列数据的最小值
min_element=min(old_data(:,1:n),[],'omitnan');   %求出第1列到第N-1列的每列数据的最小值
for i=1:n
    if max_element(:,i) ~= 0
        new_data(:,i)=(old_data(:,i)-min_element(i))/(max_element(i)-min_element(i));
    end
end

