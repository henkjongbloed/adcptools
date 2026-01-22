function plotp_compare(pvec, select_names, names, U)

% Now, pvec is a cell array of 1 x np_vec (compare np_vec vectors)


% figure;
% np = size(pvec{1},1)/mesh.ncells; % number of parameters per cell
n = length(select_names);
if numel(pvec)>1
    t = tiledlayout(n, size(pvec,2), TileSpacing="tight", Padding="tight", TileIndexing="columnmajor");
else
    t = tiledlayout('flow', TileSpacing="tight", Padding="tight");
end
t.XLabel.String = 'y [m]';
t.YLabel.String = 'z [m]';
t.XLabel.Interpreter = 'latex';
t.YLabel.Interpreter = 'latex';

for col = 1:size(pvec,2) %hardcoded
    p = p2pars(pvec{col}, U.mesh_mean.ncells);%        reshape(pvec, [np, mesh.ncells])';
    for i = 1:n
        nexttile;
        var = p(:, find(strcmp(select_names{i}, names)));
        
        loc_tit = strrep(select_names{i}, 'sig', '\sigma');
        loc_tit = strrep(loc_tit, '^1', '');
        loc_tit = strrep(loc_tit, 'u0', 'u_0');
        loc_tit = strrep(loc_tit, 'v0', 'v_0');
        loc_tit = strrep(loc_tit, 'w0', 'w_0');
        loc_tit = strrep(loc_tit, 'd', '\partial ');
        U.mesh_mean.plot(var)
        %     loc_tit = str,old,new)
        title(['$', loc_tit, '$'], 'interpreter', 'latex')
        
%         set(c,'TickLabelInterpreter','latex')
        if ~contains(loc_tit, 'phi')
            amax = max(abs(var(:,1)), [], 'omitnan') + 1e-5;
            caxis([-amax, amax])
            c=colorbar;
            colormap(gca, U.velmap)
            if ~contains(loc_tit, '\partial')
                ylabel(c, '$m/s$','Rotation',270, 'interpreter', 'latex');
            elseif contains(loc_tit, '\sigma')
                 ylabel(c, '$m/s$','Rotation',270, 'interpreter', 'latex');
            else
                 ylabel(c, '$m/s^2$','Rotation',270, 'interpreter', 'latex');
            end

        else
            caxis([-180, 180])
            temp = get(gca,'Children');
            temp(2).CData =  temp(2).CData*180/pi;
            c=colorbar;
            colormap(gca, U.phimap)
            ylabel(c, 'deg','Rotation',270, 'interpreter', 'latex');

        end
        pos = get(c,'Position');
        if i==1
            pos1 = pos;
        end
%         disp(pos)
        c.Label.Position(1) = pos1(1)+.5/col; % to change its position
%         c.Label.Position(2) = c.Label.Position(2) + .2; % to change its position

        disp(c.Label.Position)
        c.Label.HorizontalAlignment = 'center'; % to change its position
        c.TickLabelInterpreter = 'latex';
%         c.Label.Rotation = 270; % to rotate the text
        axis tight
        %     hAxes.TickLabelInterpreter = 'latex';
        %     title(sprintf('%s, %s, %s', names{i}, names{i+np(1)}, names{i + np(1) + np(2)}))
        set(gca, 'XDir','reverse') % Very important
%         xlabel('y [m]', 'interpreter', 'latex')
%         ylabel('z [m]', 'interpreter', 'latex')
        set(gca,'XTick',[])
        set(gca,'YTick',[])
    end
end
% Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
% %[Left Bottom Right Top] spacing
% NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
% set(gca, 'Position', NewPos);
end

