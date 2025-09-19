function [xns, idx, xnss, xn] = scale_sort(x)



xn = x./sum(abs(x)); % x normalized

[xns, idx] = sort(xn, 'descend', 'ComparisonMethod', 'abs'); % x sorted and with original indices 

%xns = xns.*sign(xn);
xnss = cumsum(xns); % cumulative sum of sorted normalized array

end