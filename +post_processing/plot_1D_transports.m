function plot_1D_transports(didx, bidx, SF, ri, D, bname, UU)

%% CSA f(t)

% For these figs, I use a time resulution of 4*26 instead of 26.

dim = 1;
for di = didx
    figure
    tl = tiledlayout(3,2); % One fig per dataset % two cols: SF and UU, 
    for bi = bidx
        sf0 = sum([SF{di}{bi}{ri,1}{1,:,3}]);
        nexttile()
        hold on
        plot(D{di}{bi}.X.t, sf0.*ones(size(D{di}{bi}.X.t)))
        for cidx = 1:8
            plot(D{di}{bi}.X.t, SF{di}{bi}{ri,1}{dim+1,cidx,3})
        end
        title(["$us$: ", bname{di}{bi}], 'Interpreter','latex')
        legend({"f0", SF{di}{bi}{ri,1}{dim+1,:,5}})
        hold off
        nexttile()
        uu0 = sum([UU{di}{bi}{ri,1}{1,:,3}]);
        hold on
        plot(D{di}{bi}.X.t, uu0.*ones(size(D{di}{bi}.X.t)))
        for cidx = 1:2:8
%              plot(sign(UU{di}{bi}{ri,1}{dim+1,cidx,3}).*sqrt(abs(UU{di}{bi}{ri,1}{dim+1,cidx,3}))) %plot(squeeze(UU{di}{bi}{ri,1}{4,cidx,3}), D{di}{bi}.X.sig)
                      plot(D{di}{bi}.X.t, 2*UU{di}{bi}{ri,1}{dim+1,cidx,3}) %plot(squeeze(UU{di}{bi}{ri,1}{4,cidx,3}), D{di}{bi}.X.sig)

        end
        hold off
        legend({"uu0", UU{di}{bi}{ri,1}{dim+1,1:2:8,5}} )
        title(["$u^2$: ", bname{di}{bi}]', 'Interpreter','latex')
    end
end
sgtitle("Salt Flux and $u^2$")

%% Plot f(y)

dim = 2;
for di = didx
    figure
    tl = tiledlayout(3,2); % One fig per dataset % two cols: SF and UU, 
    for bi = bidx
        sf0 = sum([SF{di}{bi}{ri,1}{1,:,3}]);
        nexttile()
        hold on
        plot(D{di}{bi}.X.y, sf0.*ones(size(D{di}{bi}.X.y)))
        for cidx = 1:8
            plot(D{di}{bi}.X.y, SF{di}{bi}{ri,1}{dim+1,cidx,3})
        end
        title(["$us$: ", bname{di}{bi}], 'Interpreter','latex')
        legend({"f0", SF{di}{bi}{ri,1}{dim+1,:,5}})
        hold off
        nexttile()
        uu0 = sum([UU{di}{bi}{ri,1}{1,:,3}]);
        hold on
        plot(D{di}{bi}.X.y, uu0.*ones(size(D{di}{bi}.X.y)))
        for cidx = 1:2:8
%              plot(sign(UU{di}{bi}{ri,1}{dim+1,cidx,3}).*sqrt(abs(UU{di}{bi}{ri,1}{dim+1,cidx,3}))) %plot(squeeze(UU{di}{bi}{ri,1}{4,cidx,3}), D{di}{bi}.X.sig)
                      plot(D{di}{bi}.X.y, 2*UU{di}{bi}{ri,1}{dim+1,cidx,3}) %plot(squeeze(UU{di}{bi}{ri,1}{4,cidx,3}), D{di}{bi}.X.sig)

        end
        hold off
        legend({"uu0", UU{di}{bi}{ri,1}{dim+1,1:2:8,5}} )
        title(["$u^2$: ", bname{di}{bi}]', 'Interpreter','latex')
    end
end
sgtitle("Salt Flux and $u^2$")



%% Plot f(sig)
dim = 3;
for di = didx
    figure
    tl = tiledlayout(3,2); % One fig per dataset % two cols: SF and UU, 
    for bi = bidx
        sf0 = sum([SF{di}{bi}{ri,1}{1,:,3}]);
        nexttile()
        hold on
        plot(sf0.*ones(size(D{di}{bi}.X.sig)), D{di}{bi}.X.sig)
        for cidx = 1:8
            plot(squeeze(SF{di}{bi}{ri,1}{dim+1,cidx,3}), D{di}{bi}.X.sig)
        end
        title(["$us$: ", bname{di}{bi}], 'Interpreter','latex')
        legend({"f0", SF{di}{bi}{ri,1}{dim+1,:,5}})
        hold off
        nexttile()
        uu0 = sum([UU{di}{bi}{ri,1}{1,:,3}]);
        hold on
        plot(uu0.*ones(size(D{di}{bi}.X.sig)), D{di}{bi}.X.sig)
        for cidx = 1:2:8
%              plot(sign(UU{di}{bi}{ri,1}{dim+1,cidx,3}).*sqrt(abs(UU{di}{bi}{ri,1}{dim+1,cidx,3}))) %plot(squeeze(UU{di}{bi}{ri,1}{4,cidx,3}), D{di}{bi}.X.sig)
                      plot(2*squeeze(UU{di}{bi}{ri,1}{dim+1,cidx,3}), D{di}{bi}.X.sig) %plot(squeeze(UU{di}{bi}{ri,1}{4,cidx,3}), D{di}{bi}.X.sig)

        end
        hold off
        legend({"uu0", UU{di}{bi}{ri,1}{dim+1,1:2:8,5}} )
        title(["$u^2$: ", bname{di}{bi}]', 'Interpreter','latex')
    end
end
sgtitle("Salt Flux and $u^2$")
end