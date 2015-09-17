function [ Days ] = Days_each_month( yr_uniqu )

for i=1:length(yr_uniqu)
    Days(i,:) = eomday(yr_uniqu(i), 1:12);
end


end

