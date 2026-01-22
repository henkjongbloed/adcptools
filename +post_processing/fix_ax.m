function ax = fix_ax(ax, opts)
%ax = nexttile;
set(ax, 'FontName',  "BookAntiqua")
set(ax, 'FontSize',  10)
axis tight
ax.TickLabelInterpreter = 'latex';
if strcmp(opts.ylim, "sym")
    YL = get(gca, 'YLim');
    maxlim = max(abs(YL));
    set(gca, 'YLim', [-maxlim maxlim]);
elseif strcmp(opts.ylim, "pos")
    YL = get(gca, 'YLim');
    maxlim = max(abs(YL));
    set(gca, 'YLim', [0 maxlim]);
end
grid on
% grid minor
if isfield(opts, "alternating_y")
    if opts.alternating_y
        c= get_c;
        nt = numel([ax.YTickLabel]);
        for i = 1:nt
            ax.YTickLabel{i} = sprintf("\\color[rgb]{%f,%f,%f}%s", ...
                c{mod(i, numel(c)) + 1}, ax.YTickLabel{i});
            disp(c{mod(i, numel(c)) + 1})
        end
    end
end
if isfield(opts, "dashedleft")
    if opts.dashedleft
        xl = xlim(); % Find out x location of the y axis.
        % Cover up existing axis with a white line.

        line([xl(1), xl(1)], ylim, 'color', 'w', 'LineWidth', 2);
        % Draw a dashed line on top of the white line.
        line([xl(1), xl(1)], ylim, 'color', 'k', 'LineStyle', '--');
    end
end
if isfield(opts, "dashedright")
    if opts.dashedright
        xl = xlim(); % Find out x location of the y axis.
        % Cover up existing axis with a white line.

        line([xl(2), xl(2)], ylim, 'color', 'w', 'LineWidth', 2);
        % Draw a dashed line on top of the white line.
        line([xl(2), xl(2)], ylim, 'color', 'k', 'LineStyle', '--','LineWidth', 2);
    end
end
if isfield(opts, "XTick")
    ax.XTick = [];
    ax.XTick = opts.XTick;
end
