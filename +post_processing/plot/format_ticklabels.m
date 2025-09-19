function fmtc = format_ticklabels(labels)

pre = "$\mathbf{||";
after = "||}$";

for idx = 1:numel(labels)
    fmtc(idx) = strjoin([pre, labels(idx), after]);
end

end

%    ["s_0", "\bar{s}_t^t", "\bar{s}_y^y", "\bar{s}_\sigma^\sigma",...
%      "\widehat{s}_{y\sigma}^{y\sigma}", "\underline{s}_{t\sigma}^{t\sigma}", "[s]_{ty}^{ty}", "s_{t y\sigma}", "s"]);
