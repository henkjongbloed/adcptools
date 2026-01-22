function plotPars(tak, U)
figure;
npars = U{1}.T.velocity_model.npars  ;
idx_0 = [1, npars(1) + 1, sum(npars(1:2)) + 1];
idx_ampl = [(idx_0(1)+1):2:(idx_0(2)-2), (idx_0(2)+1):2:(idx_0(3)-2), (idx_0(3)+1):2:sum(npars)];
idx_phase = idx_ampl + 1;
for j = 1:length(tak)
    for i = 1:size(U{j}.tid_pars,2)
        subplot(size(U{j}.tid_pars,2), length(tak), sub2ind([length(tak), size(U{j}.tid_pars,2)], j,i))
        U{j}.mesh_mean.plot_vec(U{j}.tid_pars(:,i))
        amax = max(abs(U{j}.tid_pars(:,i)), [], 'omitnan');
%         title(tit{i,j})
%         caxis([-1.5, 1.5])
        if any(i == idx_phase)
        %caxis([min(min(U{j}.tid_pars(:,d))), max(max(U{j}.tid_pars(:,d)))])
            colormap(gca, U{j}.phimap)
        else
            colormap(gca, U{j}.velmap)
%             caxis([-pi, pi])
        end
        caxis([-amax, amax])
        colorbar;
    end
end
end