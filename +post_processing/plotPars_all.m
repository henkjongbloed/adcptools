function plotPars_all(U, var)

% Ugly function hard coded for M2 and M4 (IAHR Abstract)
names = U.T.velocity_model.names;
ti = 'Case0';
p = U.pars0; % By default, unregularized parameters are plotted.


if ~isscalar(var)
    hold_stat=get(gca,'NextPlot');
    hold on
    for ce=1:numel(var)
        plotPars_all(U, var(ce))
    end
    set(gca,'NextPlot',hold_stat);
    return
end


if strcmp(var, 'Old')
    p = U.pars0;
    ti = 'Unregularized';
elseif strcmp(var, 'New')
    p = U.pars;
    ti = 'Regularized';
elseif strcmp(var, 'Difference')
    p = U.pars - U.pars0;
    ti = 'Difference';
elseif strcmp(var, 'MIOld')
    p = U.I_loc0;
    ti = 'Morans I-measure: Old';
elseif strcmp(var, 'MINEw')
    p = U.I_loc;
    ti = 'Morans I-measure: New';
end



figure('units','normalized','outerposition',[0 0 1 1])

t = tiledlayout(6,10);
sgtitle([U.BN, ' : ', ti])
for i = 1:prod(t.GridSize)
    nexttile(i);
    var = p(:,i);
    U.mesh_mean.plot(var)
    title(names{i})
    if contains(names{i}, 'phi') && contains(names{i}, '0')
        colormap(gca, U.phimap)
        caxis([-pi, pi])
    else
        colormap(gca, U.velmap)
        amax = max(abs(var(:,1)), [], 'omitnan') + 1e-3;
        caxis([-amax, amax])
    end
    colorbar;
    axis tight

    %     title(sprintf('%s, %s, %s', names{i}, names{i+np(1)}, names{i + np(1) + np(2)}))
    set(gca, 'XDir','reverse') % Very important
end
end
% if nargin > 2
%     if varargin{2}==1
%         [p, names] = ab2Aphi(p, names);
%     end
% end

% npars = U.T.velocity_model.npars  ;


% np = U.T.velocity_model.npars  ;


