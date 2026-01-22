function plotDA_Flux(tak, U, F)
% dims = {' us', ' vs'};

figure;
dim = 1;
for j = 1:length(tak)
    for d = dim
        subplot(length(tak), length(dim), sub2ind([length(tak), length(dim)], j,d))
        hold on
        for i = 1:5
            plot(U{j}.mesh_mean.n_middle, F{j}.Loc{i,d}, 'LineWidth', 2)
            grid on           
        end
        title(strcat('Transport: ', U{j}.BN), 'FontSize',12)
        axis tight
        if j==length(tak)
            xlabel('y [m]')
        end
        ylabel('psu m^2/s')
        set(gca, 'XDir','reverse')
    end
end
legend('Subtidal adv.','Tidal-sal. c-c','Stokes drift','Tidal sloshing','Vert. shear disp.', 'Location','bestoutside')

% print('lat', '-dpng')
end