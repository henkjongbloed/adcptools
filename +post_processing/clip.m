function [xclip, idxclip] = clip(x, xlim)
% input: Column vector x and limits of cutoff
% output: clipped vector and original indices that are kept

%Assumes monotonically increasing values in x.

idx = 1:numel(x);

xclip = x;
xclip(xclip < xlim(1)) = nan;
xclip(xclip > xlim(2)) = nan;

idxclip = idx(~isnan(xclip));

xclip(isnan(xclip)) = [];
end