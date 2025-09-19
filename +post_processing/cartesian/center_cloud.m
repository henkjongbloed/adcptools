function xnew = center_cloud(xold, x0)
% Center n-dimensional vector [x1, x2, ...] having column mean [xbar1,
% xbar2, ...] to predefined point x0 = [x01, x02, ...].

xbar = mean(xold,1, 'omitnan');

xnew = repmat(x0 - xbar, size(xold,1), 1) + xold;