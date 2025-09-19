function matlab_idx = matlab_idx(idx)
addone = @(c) c + 1;
if isa(idx, "double")
    matlab_idx = idx + 1;
elseif isa(idx, "cell")
    matlab_idx = cellfun(addone, idx, 'UniformOutput', false);
end

end