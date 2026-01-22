function plotTimeUS(tak, plotT, F, U)



figure;
for j = 1:length(tak)
    hours = 24*(F{j}.time - F{j}.time(1));
    for hr = 1:length(plotT) %numel(F{tak(j)}.mesh)
        subplot(numel(plotT), length(tak), sub2ind([length(tak),numel(plotT)], j, hr))
        
        vm = 0*mean(F{j}.u{plotT(hr)}(:,2), 'omitnan');
        wm = 0*mean(F{j}.u{plotT(hr)}(:,3), 'omitnan');
        csa = [(F{1, j}.uc{1, 1} + F{1, j}.uc{3, 1}(hr)), ...
            (F{1, j}.vc{1, 1} + F{1, j}.vc{3, 1}(hr)), ...
            mean(F{j}.u{plotT(hr)}(:,3), 'omitnan')];
        F{j}.mesh(plotT(hr)).plot(F{j}.u{plotT(hr)} - csa) % Look at this
        hold on
%         quiver(0, F{j}.mesh(plotT(hr)).water_level, vm, wm, "LineWidth", 3)
%         caxis([-1.5,1.5])
        hold off
        axis tight
        xlabel('n [m]')
        ylabel('z [m]')
        colormap(gca,  U{j}.velmap)
        colorbar;
        caxis([-1.5, 1.5])
        title(strcat('Flow at t = ', num2str(hours(plotT(hr))), 'h at ', U{j}.BN))
        set(gca, 'XDir','reverse') % Very important
    end
end


figure;
% suptitle('Salt')
for j = 1:length(tak)
    for hr = 1:length(plotT) %numel(F{tak(j)}.mesh)
        subplot(numel(plotT), length(tak), sub2ind([length(tak),numel(plotT)], j, hr))
        F{j}.mesh(plotT(hr)).plot(F{j}.s{plotT(hr)})
        caxis([0,max(F{j}.s{plotT(hr)}, [], 'all')]) % Conduct previous for loop not vary colorbar at different times.
        colormap(gca, U{j}.salmap)
        colorbar;
        axis tight
        title(strcat('Salt at t = ', num2str(hours(plotT(hr))), 'h at ', U{j}.BN))
        set(gca, 'XDir','reverse') % Very important
    end
end
end