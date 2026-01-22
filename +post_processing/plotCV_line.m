cvc = load("cvcoarse.mat").CV;
cvm = load("cvmod.mat").CV;
cvf = load("cvfine.mat").CV;

nreg = 30;
lc = 5;
ls = logspace(-4,3,nreg-1)';

m = makefigure(18, 6);
subplot(1,2,1)
semilogx([0;ls], [cvc{1:end,1}]/cvc{1,1})
hold on
semilogx([0;ls], [cvm{1:end,1}]/cvm{1,1})
ylim([0, max([cvc{1:end,1}]/cvc{1,1})])
legend(["Coarse", "Moderate"])
xlabel("$\lambda_s$", 'interpreter', 'latex', 'FontSize', 10);
ylabel("Scaled generalization error")
title("Coarse / moderate mesh")
xticks(ls(1:4:end))
grid on
axis tight
subplot(1,2,2)
semilogx([0;ls], [cvf{1:end,1}]/cvf{1,1})
%legend(["Fine"])
xlabel("$\lambda_s$", 'interpreter', 'latex', 'FontSize', 10);
title("Fine mesh")
axis tight
grid on
xticks(ls(1:4:end))
