function plot_1D_tidalTransport(didx, bidx, SF, ri, D, bname, UU)

%% CSA f(t)

% For these figs, I use a time resulution of 4*26 instead of 26.

dim = 1;
%makefigure(20,10)

tl = tiledlayout(2,3, "TileSpacing","compact", Padding="compact"); % One fig per dataset % two cols: SF and UU,
% tl= ["u_0", "\bar{u}_t^t", "\bar{u}_y^y", "\bar{u}_\sigma^\sigma",...
%                         "\widehat{u}_{y\sigma}", "\underline{u}_{t\sigma}", "[u]_{ty}", "u_{t y\sigma}"];
for di = didx
    for bi = bidx
        treal{di} = timeofday(datetime(D{di}{bidx(1)}.X.t, 'ConvertFrom', 'datenum'));
        tmid{di} = timeofday(datetime((D{di}{bi}.X.t(end) + D{di}{bi}.X.t(1))/2, 'ConvertFrom', 'datenum'));
        res_idx = [1, 2];
        sf0_norms = reshape([SF{di}{bi}{ri,1}{:,:,4}], [8,8]);
        resolved(di, bi) = sum(sf0_norms(1, res_idx)) + sum(sf0_norms(2, res_idx));
        unresolved(di, bi) = sum(sf0_norms(1, setdiff(1:8, res_idx))) +sum( sf0_norms(2, setdiff(1:8, res_idx)));
        perc(di, bi) = 100*resolved(di, bi)./(resolved(di, bi) + unresolved(di,bi));


        ax=nexttile();
        fontname(ax, "Book Antiqua")

        sf0 = [SF{di}{bi}{ri,1}{1,:,3}]'; %vector of numbers
        ax.FontSize = 10;
        sf0_unres = sum(sf0(setdiff(1:8, res_idx)));
        sf_tot = [sf0(res_idx); sf0_unres];

        %Norms of terms steady

        hold on
        grid on
        b=bar(tmid{di}, sum(sf0), 'stacked');
        b(1).BarWidth = 0.1;

        for cidx = 1:8 % resolved
            plot(treal{di}, SF{di}{bi}{ri,1}{dim+1,cidx,3})
        end

        %         res{di}{bi} = res{di}{bi} + sum(sf0_norms(res_idx));
        %         unres{di}{bi} = res{di}{bi} + sum(sf0_norms(setdiff(1:8, res_idx)));

        %unresolved
        %             plot(treal{di}, sum([SF{di}{bi}{ri,1}{dim+1,setdiff(1:8, res_idx),3}], 2))

        axis tight
        title([bname{di}{bi}], 'Interpreter','latex')
        if di ==2 && bi==3
            leg = legend(["$\overline{f}$", ...
                "$u_0s_1$", "$u_1s_0$", "$\overline{u_2s_6}^t$", "$\overline{u_6s_2}^t$",...
                "$\overline{u_3s_5}^t$", "$\overline{u_5s_3}^t$", "$\overline{u_4s_7}^t$", "$\overline{u_7s_4}^t$"], "Interpreter","latex", Location = "northwest", FontSize=10);
            fontsize(leg, 10, 'points')
        end
        if bi==1
            ylabel("psu m/s", 'Interpreter','latex')
        end
        %             xlabel("t [h]", 'Interpreter','latex')


    end
end
end