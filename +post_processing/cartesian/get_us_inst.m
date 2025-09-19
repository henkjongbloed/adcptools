function [U, S, D] = get_us_inst(flow, sal, evres)
% Same as get_us.m but for retrieving data for results ssec 1 in paper.




lim = {[0, flow.solver.model.periods(1,1)/(3600*24)];...
    [min(flow.solver.mesh.n_left)+.5, max(flow.solver.mesh.n_right)-.5];...
    [0,1]};
% Flow
[U,X,H] = get_var(flow, evres, lim);

% Salinity: Evaluate sigma from .1 up to .9
lim{3} = [.0, .9];
[S,~, ~] = get_var(sal, evres, lim);

D = Decomposition(X=X, H=H);

end