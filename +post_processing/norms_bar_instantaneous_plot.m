function NF = norms_bar_instantaneous_plot(D, didx, bidx, lf, NF)
[n, a] = D{didx(1)}{bidx(1)}.get_names("u");


tit = ["Flow velocity $u(t,y,\sigma)$", "Salinity $s(t,y,\sigma)$",  "Flux $f(t,y,\sigma) = (us)(t,y,\sigma)$"];
%            xticklabels(["$u_0$", "$\bar{u}_t^t$", "$\bar{u}_y^y$", "$\bar{u}_\sigma^\sigma$",...
%                "$\widehat{u}_{y\sigma}^{y\sigma}$", "$\underline{u}_{t\sigma}^{t\sigma}$", "$[u]_{ty}^{ty}$", "$u_{t y\sigma}$"]);
%            yticklabels(["$s_0$", "$\bar{s}_t^t$", "$\bar{s}_y^y$", "$\bar{s}_\sigma^\sigma$",...
%               "$\widehat{s}_{y\sigma}^{y\sigma}$", "$\underline{s}_{t\sigma}^{t\sigma}$", "$[s]_{ty}^{ty}$", "$s_{t y\sigma}$"]);
leg = {["HC", "OMS", "OMN"], ["NM", "OM", "NWW"]};
residx_ = 0:7; % change this to 0:7 for the full flow field
residx = residx_ + 1;
sgtit = ["Junction S: Inst. component norms", "Junction N: Inst. component norms"];
for di = didx
    [fig, TL] = prep_fig([10,16], [3,1]);
%     TL = tiledlayout(2,1,"TileSpacing",'compact', 'Padding','compact');
    TL.TileSpacing = "compact";
    TL.Padding = 'compact';
    for c = 1 % only orthogonal component is considered.
        % per dataset, make one aggregate plot of magnitude of components.
        for ri = 1:size(lf, 1)
            %[fig, TL] = prep_fig([20, 20], [numel(l), 4]);
%             figure
            for vi = 1:3 %flow, salt, flux
                ax = nexttile;
                dat = zeros(numel(residx)+1, numel(bidx));
                for bi = bidx
                    [NF{di}{bi}{vi}{ri, c}(:,2), NF{di}{bi}{vi}{ri, c}(:,3), NF{di}{bi}{vi}{ri, c}(:,4), NF{di}{bi}{vi}{ri, c}(:,5)] = ...
                        scale_sort(NF{di}{bi}{vi}{ri, c}(:,1));
                    % original, normalized_sorted , idx, normalized_sorted_cumsum, normalized
                    dat(1:numel(residx),bi) = sqrt(NF{di}{bi}{vi}{ri, c}(residx,1));
                    dat(numel(residx)+1,bi) = sqrt(sum(NF{di}{bi}{vi}{ri, c}(residx,1)));
                    %/
                    %                     datsf(:,bi) = SF0{di}{bi}{ri};
                end
                %                 nexttile;
                bar(residx(1):(numel(residx)+1), dat)
                if di == 1 % Junction sOUTH
                  perc = 100*sum(dat([1, 2, 3, 7], :).^2, 1)./dat(numel(residx)+1,:).^2
                else
                  perc = 100*sum(dat([1, 2, 4, 6], :).^2, 1)./dat(numel(residx)+1,:).^2
                end
                grid on
                if vi==1
                    tl= ["u_0", "\bar{u}_t^t", "\bar{u}_y^y", "\bar{u}_\sigma^\sigma",...
                        "\widehat{u}_{y\sigma}", "\underline{u}_{t\sigma}", "[u]_{ty}", "u_{t y\sigma}"];
                     tl(end+1) = "u";
%                     tl(end+1) = "\widehat{u}";
                    xticklabels(format_ticklabels(tl([residx, 9])));
                    ylabel("m/s", "FontName","Book Antiqua")
                    legend(leg{di})
                elseif vi==2
                    tl= ["s_0", "\bar{s}_t^t", "\bar{s}_y^y", "\bar{s}_\sigma^\sigma",...
                        "\widehat{s}_{y\sigma}", "\underline{s}_{t\sigma}", "[s]_{ty}", "s_{t y\sigma}"];
                     tl(end+1) = "s";
%                     tl(end+1) = "\widehat{s}";
                    
                    xticklabels(format_ticklabels(tl([residx, 9])));
                    ylabel("psu", "FontName","Book Antiqua")
                else
                    tl= ["f_0", "\bar{f}_t^t", "\bar{f}_y^y", "\bar{f}_\sigma^\sigma",...
                        "\widehat{f}_{y\sigma}", "\underline{f}_{t\sigma}", "[f]_{ty}", "f_{t y\sigma}"];
                     tl(end+1) = "f";
%                     tl(end+1) = "\widehat{s}";
                    
                    xticklabels(format_ticklabels(tl([residx, 9])));
                    ylabel("psu m/s", "FontName","Book Antiqua")
                end
                ax.TickLabelInterpreter = 'latex';
                %                 ax.LabelInterpreter = 'latex';
                ax.FontSize = 12;
                %                 set(ax, 'FontWeight', 'bold')

                %            yticklabels(["$s_0$", "$\bar{s}_t^t$", "$\bar{s}_y^y$", "$\bar{s}_\sigma^\sigma$",...
                %               "$\widehat{s}_{y\sigma}^{y\sigma}$", "$\underline{s}_{t\sigma}^{t\sigma}$", "$[s]_{ty}^{ty}$", "$s_{t y\sigma}$"]);

                %                     xticklabels(a(NF{di}{bi}{vi}{ri, c}(:,3)))
                title([tit(vi)], 'FontName','Book Antiqua', 'Interpreter','latex')
            end
            %             nexttile;
            %             bar(datsf)
            %             %             xticklabels(a(NF{1,bi, di, ri}(:,3)))
            %             title("Salt Flux")
        end

        sgtitle(sgtit(di), 'FontName', "Book Antiqua");
    end
    fontname(TL, "Book Antiqua");
end
end