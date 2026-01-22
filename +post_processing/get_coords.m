function X = get_coords(lim, res)
t = linspace(lim{1,1}(1), lim{1,1}(2), res(1)); % in days
n = linspace(lim{2,1}(1), lim{2,1}(2), res(2));
sig = linspace(lim{3,1}(1), lim{3,1}(2), res(3));

[T, N, Sig] = ndgrid(t, n, sig);

X.t = t; X.y = n; X.sig = sig;
X.T = T; X.Y = N; X.Sig = Sig;
end