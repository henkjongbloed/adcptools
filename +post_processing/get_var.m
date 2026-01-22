function evvar = get_var(var, X, reg_idx)


V = var.evaluate(reg_idx, T = X.T(:), N = X.Y(:), Sig = X.Sig(:), extrapolate=true); %Nvalues x Ndim x Nreg
evvar = cell([numel(reg_idx) ,size(V,2)]); %Nreg x Ndim
res = size(X.T);
for d = 1:size(V,2)
    for r = 1:size(V,3)
        evvar{r,d} = reshape(squeeze(V(:,d,r)), res);
    end
end


end