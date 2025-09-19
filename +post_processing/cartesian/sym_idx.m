function sym_idx = sym_idx(idx)
% input: cell of indices
% output: cell of indices and their swapped counterparts
sym_idx = cell([size(idx,1), 2*size(idx, 2)]);
for i = 1:size(idx,1)
    for j = 1:size(idx, 2)
        [sym_idx{i, 2*j-1}] = idx{i,j};
        [sym_idx{i, 2*j}] = fliplr(idx{i, j});
    end
end
end

