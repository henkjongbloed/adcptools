function [fig, TL] = prep_fig(dim, tldim)
    fig  = makefigure(dim(1), dim(2));
%     fig = figure;
    TL = tiledlayout(tldim(1), tldim(2));
%     set(fig, 'DefaultAxesFontSize', 10);
%     set(fig, 'DefaultAxesFontSize', 10);
%     font = 'BookAntiqua';
%     %set('DefaultTextFontName', font);
%    % set('DefaultAxesFontName', font);

end
