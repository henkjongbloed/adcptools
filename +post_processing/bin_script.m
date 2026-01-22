% bin

%% Closer investigation of H - swap y axis? - Yes
% 
% % Xs orientations
% 
% for di = didx
%  figure;
%     for bi = bidx
% %         nexttile;
% %         flow{di}{bi}.solver.xs.plot();
% %         [DH{di}{bi}, AH{di}{bi}] = D{di}{bi}.decompose_function(D{di}{bi}.H);
% %         D{di}{bi}.plot_components( AH{di}{bi}, "salmap")
% %         sgtitle(["H - avg", bname{di}{bi}])
% % 
% %         D{di}{bi}.plot_components( DH{di}{bi}, "salmap")
% %         sgtitle(["H - orth", bname{di}{bi}])
%                 nexttile;
%                 surf(squeeze(D{di}{bi}.H(:,:,1)))
%                 xlabel("y")
%                 ylabel("t")
%                 colorbar;
%                 title(bname{di}{bi})
%     end
% end
% 
% 
% for di = didx
%     figure;
%     for bi = bidx
%         hold on
%         flow{di}{bi}.solver.xs.plot();
%     end
% end
%%

% Inspect H
% for di = didx
%     figure;
%     for bi = bidx
%         flow{di}{bi}.solver.xs.plot();
%         [DH{di}{bi}, AH{di}{bi}] = D{di}{bi}.decompose_function(D{di}{bi}.H);
%         D{di}{bi}.plot_components( AH{di}{bi}, "salmap")
%         sgtitle(["H - avg", bname{di}{bi}])
% 
%         D{di}{bi}.plot_components( DH{di}{bi}, "salmap")
%         sgtitle(["H - orth", bname{di}{bi}])
%         %         nexttile;
%         %         surf(squeeze(D{di}{bi}.H(:,:,1)))
%         %         xlabel("y")
%         %         ylabel("t")
%         %         colorbar;
%         %         title(tak_names{di}{bi})
%         % disp(mean(D{di}{bi}.H,"all"))
%     end
% end



% D{di}{bi}.plot_components(cellfun(@times, DU{di}{bi}, DS{di}{bi}, 'UniformOutput',false), 'salmap')
% sgtitle(["Transport: ", tak_name{di}{bi}])

%
% figure
% tak=1;
% for t = 1:evres(1)
%     subplot(3,2,2*ti-1)
%     hold on
%     for t = 1:evres(1)
%         plot(D{ti}.X.y, squeeze(U{ti}{1}(t, :, 3)), 'b')
%         plot(D{ti}.X.y, squeeze(U{ti}{1}(t, :, 17)), 'r')
%     end
%     hold off
% end
%
%
%
%
% M2T = flow.solver.model.periods(1,1)/(3600*24);
% visres = [12, evres(2), 3];
% t = linspace(0, M2T, visres(1)); % in days
% n = linspace(min(flow.solver.mesh.n_left)+.5, max(flow.solver.mesh.n_right)-.5, visres(2));
% sig = linspace(.2, .8, visres(3));
%
%
%
%
% DF=D.decompose_function(U);
% D.plot_components(DF, "velmap");
% DFs=D.decompose_function(S);
% D.plot_components(DFs, "salmap");
%
% [SF, ~] = D.decompose_product(U,U);
% figure;
% bar(SF)
% [SF, ~] = D.decompose_product(U,S);
% figure;
% bar(SF)

% %% Map figures
%
% %small: dataset 1
% plotMap(dataset_name{1}, 'small')
% plotMap_xs(dataset_name{1}, 'small')
% %plotMap_xs(dataset_name{2}, 'large')
%
%
% %small: dataset 2
% figure
% plotMap(dataset_name{2}, 'small')
% plotMap_xs(dataset_name{2}, 'small')
% %plotMap_xs(dataset_name{2}, 'large')
%
%
% %large: dataset 1 and 2
% plotMap(dataset_name{1}, 'large')
% plotMap_xs(dataset_name{1}, 'large')
% plotMap_xs(dataset_name{2}, 'large')
%
%
% % Decomposition, visualization and rendering paper figures.
% %D.plot_components(DF, "velmap")
% %D.plot_components(DFs, "salmap")
%
% %flow.evaluate(1, )
%
% % [DF, AF, DFf, AFf] = D.decompose_function(u);
%
%
% %flow.plot_solution({flow_model.all_names{[1]}}, representation = "Aphi", sol_idx = 1) % u,v,w
% %sal.plot_solution({sal_model.all_names{[1]}}, representation = "Aphi", sol_idx = 1) % u,v,w
% %clim([0,15])
% %colormap(helpers.cmaps("salmap"))
%
% %Dsums = sum(cat(4, DFfs{:}), 4); %time x lateral x vertical
%
% %animate_solution(Dsum, X, nqt, Dsums);
%
% % Top-view bathymetries
% plotMap(dataset_name{2}, 'small')
% [xr, xl] = utm2ll(x,y, 31, 'wgs84');
% geoplot(xr,xl)
% geoscatter
% for di = 2
%     %[adcp, ctd, h, dataset_name{di}, ~] = import_data2(di);
%     %figure;
%     hold on
%     for ti = 1:3
%         flow{ti, di}.solver.bathy.plot_dots()
%         flow{ti, di}.solver.xs.plot()
%         %figure
%         %D{ti, di}.plot_components(DF{ti, di}, 'velmap')
%         %figure
%         %D{ti, di}.plot_components(DFs{ti, di}, 'salmap')
%
%         %D{ti, di}.plot_components(AF{ti, di}, 'velmap')
%         %figure
%         %D{ti, di}.plot_components(AFs{ti, di}, 'salmap')
%     end
% end
% hold off
% % Waterlevels
%
% %meshes: Fig 1 and 2
% for di = 1
%     %[adcp, ctd, h, dataset_name{di}, ~] = import_data2(di);
%     for ti = 1:3
%         m = makefigure(10,3);
%         flow{ti, di}.plot_mesh(gca);
%         disp("max depth")
%         min(flow{ti, di}.solver.bathy.known(3,:) )
%         disp("mean depth")
%         mean(flow{ti, di}.solver.bathy.known(3,:) )
%         disp("width ")
%         sum(abs(flow{ti, di}.solver.mesh.nw  ))
%     end
% end


%% bars: decomposition of fluxes forming the 1D tidal transport




% figure;
% vn = "t";
% for di = didx
%     for bi = bidx
%         for ri = 1:size(lf, 1)
%             [SF0t{di}{bi}{ri,1}, SF1t{di}{bi}{ri,1}] = D{di}{bi}.decompose_product_bar(vn, DF{di}{bi}{1}{ri, 1}, DF{di}{bi}{2}{ri, 1}); % Salt Flux
%             nexttile;
%             hold on
%             title(bname{di}{bi})
%             [~, idxd] = D{di}{bi}.get_bar(vn);
%             nm = ["00", "11", "22", "33", "44", "55", "66", "77",...
%                 idxd{:}];
%             for i = 1:8
%                 plot(D{di}{bi}.X.t, SF0t{di}{bi}{ri,1}{i}.*ones(size(D{di}{bi}.X.t)), LineWidth=2)
%             end
%
%
%             idxd = cellfun(@num2str, idxd, 'UniformOutput',false);
%             for i = 1:8
%                 plot(D{di}{bi}.X.t, SF1t{di}{bi}{ri,1}{i}, 'LineStyle','--', LineWidth=2, Tag=nm{8+i})
%             end
%             if di==2 && bi==3
%             legend(nm)
%             end
%         end
%     end
% end




%% bars: decomposition of fluxes forming the 1DH subtidal transport

% figure;
% vn = "y";
% for di = didx
%     for bi = bidx
%         for ri = 1:size(lf, 1)
%             [SF0y{di}{bi}{ri,1}, SF1y{di}{bi}{ri,1}] = D{di}{bi}.decompose_product_bar(vn, DF{di}{bi}{1}{ri, 1}, DF{di}{bi}{2}{ri, 1}); % Salt Flux
%             nexttile;
%             hold on
%             title(bname{di}{bi})
%             for i = 1:8
%                 plot(D{di}{bi}.X.y, SF0y{di}{bi}{ri,1}{i}.*ones(size(D{di}{bi}.X.y)), LineWidth=2)
%             end
%             [~, idxd] = D{di}{bi}.get_bar(vn);
%             idxd = cellfun(@num2str, idxd, 'UniformOutput',false);
%             for i = 1:8
%                 plot(D{di}{bi}.X.y, SF1y{di}{bi}{ri,1}{i}, 'LineStyle','--', LineWidth=2)
%             end
%             if di==2 && bi==3
%             legend(["00", "11", "22", "33", "44", "55", "66", "77",...
%                 idxd{:}])
%             end
%         end
%     end
% end



%% bars: decomposition of fluxes forming the 1DV subtidal transport
% figure;
% vn = "sig";
% for di = didx
%     for bi = bidx
%         for ri = 1:size(lf, 1)
%             [SF0s{di}{bi}{ri,1}, SF1s{di}{bi}{ri,1}] = D{di}{bi}.decompose_product_bar(vn, DF{di}{bi}{1}{ri, 1}, DF{di}{bi}{2}{ri, 1}); % Salt Flux
%             nexttile;
%             hold on
%             title(bname{di}{bi})
%             for i = 1:8
%                 plot(squeeze(SF0s{di}{bi}{ri,1}{i}).*ones(size(D{di}{bi}.X.sig)), D{di}{bi}.X.sig, LineWidth=2)
%             end
%             [~, idxd] = D{di}{bi}.get_bar(vn);
%             idxd = cellfun(@num2str, idxd, 'UniformOutput',false);
%             for i = 1:8
%                 plot(squeeze(SF1s{di}{bi}{ri,1}{i}), D{di}{bi}.X.sig, 'LineStyle','--', LineWidth=2)
%             end
%             if di==2 && bi==3
%             legend(["00", "11", "22", "33", "44", "55", "66", "77",...
%                 idxd{:}])
%             end
%         end
%     end
% end

%% test speed of decomposition procedure %full decomposition approx .2 sec

% temp = D{2}{1}.avg_op_all(DF{2}{1}{1}{1, 1}, DF{2}{1}{2}{1, 1}, {"y", "sig"});
% % D{2}{1}.plotC(cell2mat(temp))
%


