function quick_inspect_plot(didx, bidx, lf, D, DF, AF, c)

cm = ["velmap", "salmap", "fluxmap"];

for di = didx
    for bi = bidx
        for ri = 1:size(lf, 1)
            for vi = 1:2 %flow, salt, flux
%                 figure;
%                 D{di}{bi}.plot_components(AF{di}{bi}{vi}{ri, 1}, cm(vi))
%                 sgtitle(["Averaged ", tit(vi), bname{di}{bi}, "with l = ", num2str(l(ri))])

                D{di}{bi}.plot_components(DF{di}{bi}{vi}{ri, 1}, cm(vi))
%                 sgtitle(["Orthogonal ", tit(vi), bname{di}{bi}, "with l = ", num2str(l(ri))])
            end
        end
    end
end
end