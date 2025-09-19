function plotTZ(U, F, nnorm)

% Plot velocity and salinity in time and in the vertical, for a fixed
% lateral coordinate n.

br = numel(F);
m = zeros(4,1);
t0 = 24*F{1}.time(1);
for j = 1:br
%     t{j} = F{j}.time' ;
    t{j} = 24*F{j}.time' - t0;
    n = nnorm*U{j}.Bw/2;
    for i = 1:length(t{j})
        sig{j}(:,i) = linspace(.1, .9, 10);
        Z{j}(:,i) = F{j}.mesh(i).sigma_to_z(sig{j}(:,i), n2zb(n, F{j}.mesh(i)), F{j}.eta(i));
        ind = F{j}.mesh(i).index(n*ones(size(sig{j}(:,i))), sig{j}(:,i));
        for nc = 1:length(ind) % Ugly but works
            if isnan(ind(nc))
                ind(nc) = F{j}.mesh(i).index(n*ones(size(sig{j}(nc,i))), sig{j}(nc-1,i));
            end
        end
        u{j}(:,i) = F{j}.u{1,i}(ind,1);
        v{j}(:,i) = F{j}.u{1,i}(ind,2);
        w{j}(:,i) = F{j}.u{1,i}(ind,3);
        s{j}(:,i) = F{j}.s{1,i}(ind,1);

    end
    T{j} = repmat(t{j}, size(sig{j}, 1), 1);
    m(1) = max(m(1),max(abs(u{j}(:))));
    m(2) = max(m(2),max(abs(v{j}(:))));
    m(3) = max(m(3),max(abs(w{j}(:))));
    m(4) = max(m(4),max(abs(s{j}(:))));
end
% c = sin(3600*T);


figure;
varName = {'u', 'v', 'w', 's'};
for j = 1:br
    vars = {u{j}, v{j}, w{j}, s{j}};
    for d = 1:4
        ax=subplot(4,br,sub2ind([br,4], j,d));
        contourf(T{j}, Z{j}, vars{d}, 50, 'edgecolor', 'none');
        shading interp;
        hold on
        plot(t{j}, F{j}.eta, 'LineWidth', 2)

        shading interp;
        %set(c, 'edgecolor','none');
        colormap(gca, U{j}.velmap);
        caxis([-1, 1].*m(d))
        if d == 4
           colormap(gca, U{j}.salmap);
           caxis([0, 1].*m(d))
        end
        colorbar;
        title([U{j}.BN, ': ', varName{d}])
        xlabel('t [hour]')
        ylabel('z [m]')
        xticks(ax, [0:12]');
        %xticklabels({'0', [], '2', [], '4', [], '6', [], '8', [], '10', [], '12'})
        ax.XGrid = 'on';
    end
end

end