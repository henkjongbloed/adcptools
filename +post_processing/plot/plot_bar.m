function plot_bar(didx, bidx, SF, ri, dim, bname, UU)
tl = tiledlayout(3,2); % One fig per dataset % two cols: SF and UU, 
for di = didx
    for bi = bidx
        nexttile()
        hold on
        for cidx = 1:8
            plot(SF{di}{bi}{ri,1}{dim,cidx,3})
        end
        title(["u^2: ", bname{di}{bi}])
        legend(SF{di}{bi}{ri,1}{dim,:,5} )
        nexttile()
        for cidx = 1:2:8
            plot(UU{di}{bi}{ri,1}{dim,cidx,3})
        end
        legend(SF{di}{bi}{ri,1}{dim,:,5} )
        title(["u^2: ", bname{di}{bi}])
    end
end
end