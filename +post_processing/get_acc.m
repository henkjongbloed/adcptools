function [acc,X] = get_acc(flow, sol_idx, t0)

% Function that extract the acceleration term from the solved velocity.

nqt = 2; nqn = 50; nqs = 20;

secday = 3600*24;
dt = 1/(24*60); % 1 minute
t = linspace(t0-dt, t0+dt, nqt);
n = linspace(min(flow.solver.mesh.n_left)+.5, max(flow.solver.mesh.n_right)-.5, nqn);%*ones(size(t));
dn = n(2)-n(1);
sig = linspace(0, 1, nqs);%*ones(size(t));
dsig = sig(2)-sig(1);
dx = .5;
[T, N, Sig] = ndgrid(t,n,sig); Tq = T(:); Nq = N(:); Sigq = Sig(:);

U{1} = flow.evaluate(sol_idx, T = Tq, N = Nq, Sig = Sigq, dX = dx, extrapolate=true);
U{2} = flow.evaluate(sol_idx, T = Tq, N = Nq, Sig = Sigq, dX = 0, extrapolate=true);
U{3} = flow.evaluate(sol_idx, T = Tq, N = Nq, Sig = Sigq, dX = -dx, extrapolate=true);

% [T, N, Sig] = ndgrid(t,n,sig); Tq = T(:); Nq = N(:); Sigq = Sig(:);
u{1} = reshape(U{1}(:,1), [nqt, nqn, nqs]);
v{1} = reshape(U{1}(:,2), [nqt, nqn, nqs]);
w{1} = reshape(U{1}(:,3), [nqt, nqn, nqs]);

u{2}= reshape(U{2}(:,1), [nqt, nqn, nqs]);
v{2} = reshape(U{2}(:,2), [nqt, nqn, nqs]);
w{2} = reshape(U{2}(:,3), [nqt, nqn, nqs]);

u{3}= reshape(U{3}(:,1), [nqt, nqn, nqs]);
v{3} = reshape(U{3}(:,2), [nqt, nqn, nqs]);
w{3} = reshape(U{3}(:,3), [nqt, nqn, nqs]);


% Determine dw/dx
dwdx = squeeze(mean(w{1} - w{3}, 1));
dwdt = squeeze((w{2}(2,:,:) - w{2}(1,:,:))/(2*dt))/(secday);
u0 = squeeze(mean(u{2}, 1));
v0 = squeeze(mean(v{2}, 1));
w0 = squeeze(mean(w{2}, 1));

[dwdy, dwdsig] = gradient(w0, dn, dsig);

X.t = t; X.y = n; X.sig = sig;
X.T = T; X.Y = N; X.Sig = Sig;

[Hf, Wlf, Zbf] = get_H(X, flow.solver.mesh, flow.solver.bathy, flow.solver.adcp.water_level_object, 1); % Own function to derive H(t,y) from X and adcptools variables.

H = repmat(Hf, [1,1,numel(X.sig)]);
H0 = squeeze(mean(H, 1));

Zb = repmat(Zbf, [1,1,numel(X.sig)]); 
X.Z = Zbf + X.Sig.*H;

acc{1} = dwdt;
acc{2} = u0.*dwdx;
acc{3} = v0.*dwdy;
acc{4} = w0.*dwdsig./H0;
end

% fi = figure;
% filename = 'testAnimated.gif';
% ncolor = 100;
% amax = max(abs([min(Dsum, [], "all"), max(Dsum, [], "all")]));
% levels = linspace(-amax, amax, ncolor);
% [~,ha]=contourf(squeeze(X.Y(1,:,:))', squeeze(X.Z(1,:,:))', squeeze(Dsum(1,:,:))' , levels, "LineColor",'none');
% colorbar;
% title("Hartelkanaal: Flow")
% colormap(gca, flipud(brewermap(ncolor, 'RdBu')));
% clim([-amax, amax])
% ylim([min(X.Z, [], 'all'), max(X.Z, [], 'all')])
% set(gca, 'XDir','reverse') % Very important
