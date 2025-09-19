function plotPars_all_quiver(U, varargin)

% Ugly function hard coded for M2 and M4 (IAHR Abstract)
ti = 'Case0';
p = U.pars0;
if nargin>1
    if varargin{1} == 0 
        p = U.pars0;
        ti = 'CaseOld';
    elseif varargin{1} == 1 
        p = U.pars;
        ti = 'CaseNew';
    elseif varargin{1} == -1
        p = U.pars - U.pars0;
        ti = 'CaseDiff';
end
figure;
% npars = U.T.velocity_model.npars  ;
names = U.T.velocity_model.names;

np = U.T.velocity_model.npars  ;

t = tiledlayout(4,5);
sgtitle(ti)
for i = 1:prod(t.GridSize)
    nexttile(i);
    var = p(:,[i, i+np(1), i + np(1) + np(2)]);
    U.mesh_mean.plot(var)
    colormap(gca, U.velmap)
    amax = max(abs(var(:,1)), [], 'omitnan') + 1e-3;
    caxis([-amax, amax])
    colorbar;
    axis tight
    title(sprintf('%s', names{i}))
%     title(sprintf('%s, %s, %s', names{i}, names{i+np(1)}, names{i + np(1) + np(2)}))
    set(gca, 'XDir','reverse') % Very important
end

end