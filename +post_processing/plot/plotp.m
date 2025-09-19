function plotp(pvec, select_names, names, n_col, U)

% Ugly function hard coded for M2 and M4 (IAHR Abstract)
% figure('units','normalized','outerposition',[0 0 1 1])
% figure;
np = length(pvec)/U.mesh_mean.ncells; % number of parameters per cell
p = reshape(pvec, [np, U.mesh_mean.ncells])';
n = length(select_names);
n_row = ceil(n/n_col);
if n<n_col
    n_col = n;
    n_row = 1;
end
t = tiledlayout(n_row,n_col, TileSpacing="tight", Padding="tight");
% sgtitle([U.BN, ' : ', ti])
for i = 1:n
    nexttile;
    var = p(:, find(strcmp(select_names{i}, names)));
    U.mesh_mean.plot(var)
    loc_tit = strrep(select_names{i}, 'sig', '\sigma');
    loc_tit = strrep(loc_tit, '^1', '');
    loc_tit = strrep(loc_tit, 'u0', 'u_0');
    loc_tit = strrep(loc_tit, 'v0', 'v_0');
    loc_tit = strrep(loc_tit, 'w0', 'w_0');
    loc_tit = strrep(loc_tit, 'd', '\partial ');

    %     loc_tit = str,old,new)
    title(['$', loc_tit, '$'], 'interpreter', 'latex')

    if ~contains(loc_tit, 'phi')
        colormap(gca, U.velmap)
        amax = max(abs(var(:,1)), [], 'omitnan') + 1e-5;
        caxis([-amax, amax])
    else
        colormap(gca, U.phimap)
        caxis([-pi, pi])
    end
    c=colorbar;
    set(c,'TickLabelInterpreter','latex')
    %     ylabel(c, 'm/s','Rotation',270, 'interpreter', 'latex');
    axis tight
    %     hAxes.TickLabelInterpreter = 'latex';
    %     title(sprintf('%s, %s, %s', names{i}, names{i+np(1)}, names{i + np(1) + np(2)}))
    set(gca, 'XDir','reverse') % Very important
    xlabel('y [m]', 'interpreter', 'latex')
    ylabel('z [m]', 'interpreter', 'latex')
    set(gca,'XTick',[])
    set(gca,'YTick',[])
end
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);
end

