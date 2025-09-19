function residual_top_plot(didx, bidx, flow, salt,F, D, AF)
tit = ["Junction S: Residual dynamics", "Junction N: Residual dynamics"];


llimpos = [.3, .4];
llimneg = [-.4, -.3];


vlimbot = [.1, .2];
vlimsur = [.8, .9];


for di = didx
    for bi = bidx
        for vi = 1:3 %Flow, Salt, Flux
            for ri = 1
                for c = 1:size(F{di}{bi}{vi}, 2)

                    % LATERAL SHEAR VARIABLES
                    IFP{di}{bi}{vi}{ri, c} = D{di}{bi}.extract_inst_llim(F{di}{bi}{vi}{ri, c}, llimpos);
                    IFN{di}{bi}{vi}{ri, c} = D{di}{bi}.extract_inst_llim(F{di}{bi}{vi}{ri, c}, llimneg);

                    daIFP{di}{bi}{vi}{ri, c} = D{di}{bi}.da(IFP{di}{bi}{vi}{ri, c});
                    daIFN{di}{bi}{vi}{ri, c} = D{di}{bi}.da(IFN{di}{bi}{vi}{ri, c});


                    % BOTTOM AND SURFACE VARIABLES
                    IFB{di}{bi}{vi}{ri, c} = D{di}{bi}.extract_inst_vlim(F{di}{bi}{vi}{ri, c}, vlimbot); % U-Flow, Salinity, Flux instantaneous
                    IFS{di}{bi}{vi}{ri, c} = D{di}{bi}.extract_inst_vlim(F{di}{bi}{vi}{ri, c}, vlimsur); % U-Flow, Salinity, Flux instantaneous

                    % Laterally averaged bottom and surface variables.
                    laIFB{di}{bi}{vi}{ri, c} = D{di}{bi}.la(IFB{di}{bi}{vi}{ri, c});
                    % laterally averaged stratification for all components, function of (t)
                    laIFS{di}{bi}{vi}{ri, c} = D{di}{bi}.la(IFS{di}{bi}{vi}{ri, c});

                    % Tidally averaged bottom and surface variables.
                    twaIFB{di}{bi}{vi}{ri, c} = D{di}{bi}.extract_inst_vlim(AF{di}{bi}{vi}{ri, c}{5}, vlimbot); % TWA of the above
                    twaIFS{di}{bi}{vi}{ri, c} = D{di}{bi}.extract_inst_vlim(AF{di}{bi}{vi}{ri, c}{5}, vlimsur); % TWA of the above
                end
            end
        end
    end
end







for di = didx
    [fig, TL] = prep_fig([10,10], [1,1]);
    ax = plot_dynamics_residual(di, bidx, flow, salt, D, AF, twaIFB, twaIFS);
    axis tight
    title(tit(di))
    fontname(TL, "Book Antiqua")
%     set(ax, 'visible', 'off')
    ax.XTickLabel = [];
    ax.XTick = [];
    ax.YTickLabel = [];
    ax.YTick = [];
    topxlim{di} = ax.XLim;
    topylim{di} = ax.YLim;
%
%    topxlim{di} = ax.XLim;
%     topylim{di} = ax.YLim;
end

%% Plot tidal dynamics with top view

tplot = 1:4:24; % indices, so integers >0 -> approx 2 hours between each plot.
titi{1} = 4:2:14;
titi{2} = 5:2:15;
for di = didx
    [fig, TL] = prep_fig([16.4,10], [2,3]);
    %     figure;
    i = 1;
    %tpp=0;
    for ti = tplot
%         tit = datetime(D{di}{bidx(1)}.X.t(ti)  , 'ConvertFrom', 'datenum');
        ax=nexttile;
        %tpp = tpp+1;
        for bi = bidx
            [ax, maxuu(di, bi, i), maxss(di,bi,i)] = plot_dynamics_instantaneous(ax, flow, salt, D, AF, IFB, IFS, di, bi, ti);
        end
        axis square
        title(sprintf("t = %i:00 h", titi{di}(i)))
        i = i+ 1;
        fontname(TL, "Book Antiqua")
%             set(ax, 'visible', 'off')
        ax.XTickLabel = [];
        ax.XTick = [];
        ax.YTickLabel = [];
        ax.YTick = [];
        set(ax, "XLim", topxlim{di} + 50*[-1, 1])
        set(ax, "YLim", topylim{di} + 50*[-1, 1])
    end
    TL.TileSpacing = 'tight';
    TL.Padding = 'tight';
    %linkaxes([ax{:}],'x')
%     TL.XLabel.String = 't [h]';
%     %     t.YLabel.String = 'z';
%     TL.XLabel.Interpreter = 'latex';
%     TL.YLabel.Interpreter = 'latex';
    fontname(TL, "Book Antiqua")
end
end