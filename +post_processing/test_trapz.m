% test-trapz

t = 0:.001:2*pi;
f = @(t) sin(t);

F = f(t);


It = @(t, f) trapz(f(t), t);
Im = @(t, f) mean(f(t), 'all');


It(t,f)
Im(t,f)