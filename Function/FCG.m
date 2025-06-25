function [card] = FCG(gb_list_final, gb_card, GB_num, n, Q)


sumcard = 0;
for r = 1: GB_num
    gbnum = size(gb_list_final{r},1);
    temp = gbnum -  gb_card(r)/n ;
%     sumcard = sumcard + temp / gbnum;
       sumcard = sumcard + temp * Q(r) / gbnum;
end

%sumcard * Q 
card = sumcard / GB_num;
end
