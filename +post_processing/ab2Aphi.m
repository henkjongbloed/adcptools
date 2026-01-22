function [tid_pars, tid_names] = ab2Aphi(pars, names)
% Function that converts tidal constituents and gradients to A - phi
% representation
tid_pars = nan(size(pars));
tid_names = names;
for nr = 1:size(pars,3)
    ind = 1;
    while ind <= length(names) % Not vectorized still.
        if contains(names{ind}, 'M0') % Subtidal, spatially constant and varying
            tid_pars(:,ind,nr) = pars(:,ind,nr);
            %tid_names{1, ind} = strrep(names{1,ind}, 'M0A', 'M0');
            ind = ind + 1;
        elseif ~contains(names{ind}, 'd') % Tidal, spatially constant
            tid_pars(:,ind,nr) = sqrt(pars(:,ind,nr).^2 + pars(:,ind+1,nr).^2);
            tid_pars(:,ind+1,nr) = atan2(pars(:,ind+1,nr),pars(:,ind,nr));
            tid_names{1, ind} = strrep(names{1,ind}, 'a', 'A');
            tid_names{1, ind+1} = strrep(names{1,ind+1}, 'b', '\phi');
            ind = ind + 2;
            if ind <= length(names)
                 tempa = pars(:,ind,nr); tempb = pars(:,ind+1,nr);
            end
        else % Tidal, spatially varying -> Chain rule
            tid_pars(:,ind,nr) = (tempa.*pars(:,ind,nr) + tempb.*pars(:,ind+1,nr))./(sqrt(tempa.^2 + tempb.^2));
            tid_pars(:,ind + 1,nr) = (tempa.*pars(:,ind+1,nr) - tempb.*pars(:,ind,nr))./(tempa.^2 + tempb.^2);
            tid_names{1, ind} = strrep(names{1,ind}, 'a', 'A');
            tid_names{1, ind+1} = strrep(names{1,ind+1}, 'b', '\phi');
            ind = ind+2;
        end
    end
end


end