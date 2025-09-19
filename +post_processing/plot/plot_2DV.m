function [resolved, unresolved, perc] =  plot_2DV(fig, didx, bidx, SF, ri, D, bname, resolved, unresolved, perc)
dim = 3;
% makefigure(10,10)
tl = tiledlayout(2,1, "TileSpacing","compact", Padding="compact");% One fig per dataset % two cols: SF and UU,
% tl= ["u_0", "\bar{u}_t^t", "\bar{u}_y^y", "\bar{u}_\sigma^\sigma",...
%                         "\widehat{u}_{y\sigma}", "\underline{u}_{t\sigma}", "[u]_{ty}", "u_{t y\sigma}"];
for di = didx
    for bi = bidx
        % select only OMN and NWW
        sf0_norms = reshape([SF{di}{bi}{ri,1}{:,:,4}], [8,8]);
            res_idx = [1, 4];
        resolved(di, bi) = sqrt(sum(sf0_norms(1, [1,4])) + sum(sf0_norms(4, 1:2)));
        unresolved(di, bi) = sqrt(sum(sf0_norms(1, setdiff(1:8, res_idx))) + sum(sf0_norms(4, setdiff(1:8, 1:2))));
        perc(di, bi) = 100*resolved(di, bi)./(resolved(di, bi) + unresolved(di,bi));
        if bi==3

            ax=nexttile();
            fontname(ax, "Book Antiqua")
            ax.FontSize = 10;
            
            sf0 = [SF{di}{bi}{ri,1}{1,:,3}]'; %vector of numbers

            sf0_unres = sum(sf0(setdiff(1:8, res_idx)));
            sf_tot = [sf0(res_idx); sf0_unres];

            hold on

            b = barh(.5, sf_tot, 'stacked');
            b(1).BarWidth = 0.1;

            for cidx = 1:2 % resolved
                plot(squeeze(SF{di}{bi}{ri,1}{dim+1, cidx, 3}), D{di}{bi}.X.sig)
            end

            %unresolved
            plot(sum(squeeze([SF{di}{bi}{ri,1}{dim+1, setdiff(1:8, 1:2), 3}]),1), D{di}{bi}.X.sig)

            axis tight
            title(strcat("$\overline{us}^\sigma$: ", bname{di}{bi}), 'Interpreter','latex')
            if di == 2
                legend(["$u_0s_0$", "$[u_3s_3]$", "Unresolved", "$u_0s_3$", "$u_3s_0$",  "Unresolved"], "Interpreter","latex", "Location","southwest")
            xlabel("psu m/s", 'Interpreter','latex')
            end
                grid on
            ylabel("$\sigma$", 'Interpreter','latex')
            
        end
    end
end
end