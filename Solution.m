classdef Solution < handle & helpers.ArraySupport
    %Solution Class for capturing model parameters as output of
    %get_parameters

    %   Detailed explanation goes here

    properties(SetAccess = ?Solver)
        M (:,:) double % matrix M such that b = Mp

        b (:,1) double % rhs of system of eqs (b = Mp)

        p (:,:) double % model parameters (b = Mp)

        cell_idx (:,1) double

%         mesh (1,1) Mesh = SigmaZetaMesh
% 
%         regularization (1,:) regularization.Regularization
% 
%         model (1,1) DataModel = TaylorScalarModel

        solver (1,1)

        opts (1,1) SolverOptions

       % ns (:,:) double

        GOF % Goodness of fit results

    end

        %  Solver properties:
    %   adcp - VMADCP object with the adcp data
    %   mesh - Mesh object defining a mesh
    %   bathy - Defines the bathymetry
    %   xs - Defines the cross-section
    %   ensemble_filter - defines repeat transects to be processed
    %   data_model - defines the velocity model to fit

    properties
        decomp (1,1) Decomposition
    end
    properties(Dependent)
        pars %p arranged in matrix of size Nc x np.
        ns   %number of measurements per cell.
    end
    methods
        function obj = Solution(varargin)
            obj = obj@helpers.ArraySupport(varargin{:})
        end
        function pars = get.pars(obj)
            pars = obj.p2pars(obj.p);
        end
        function ns = get.ns(obj)
            ns = accumarray(obj.cell_idx, ones(size(obj.cell_idx)),[obj.solver.mesh.ncells, 1]);
        end
        function varargout = get_data(obj)
            varargout = cell(1,nargout);
            [varargout{:}] = obj.solver.model.get_data(obj.pars);
        end

        function gof = validate(obj)
            if ~isscalar(obj)
                obj.run_method('validate');
                return
            end
            gof = obj.get_residuals();
        end

        function res = evaluate(obj, sol_idx, X)
            % Wrapper for model.get_data, with arbitrary number of query
            % points
            % input: sol_idx: index of the solution to be evaluated
            % X.T : Nq x 1 vector of time values
            % X.N : Nq x 1 vector of lateral coordinate values
            % X.Sig : Nq x 1 vector of sigma values

            % evaluate(obj, sol_idx, T = ..., N = ..., Sig = ...)
            % Used for interpolation or plotting
            arguments
                obj
                sol_idx (1,1) double  = 1
                X.T (:,1) double = 0
                X.dX (1,1) double = 0 % scalar used to investigate along-channel gradients.
                X.N (:,1) double = 0
                X.Sig (:,1) double = 0
                X.extrapolate (1,1) logical = false
            end
            % if statements: numel XT, XN, XSig must be equal
            % npars must have equal elements
            nx = numel(X.T);
            np = obj.solver.model.npars(1);

            % find mesh cell query points belong to
            ev_cidx = obj.solver.mesh.index(X.N, X.Sig);

            % Query points too close to bottom or surface.

            if X.extrapolate
                ev_cidx = obj.solver.mesh.index_extrapolate(ev_cidx, X.N, X.Sig);
            else
                error("Query points outside mesh domain. " + ...
                    "Enter extrapolate = true to linearly " + ...
                    "extrapolate the solution.")
            end
            n_center = reshape(...
                obj.solver.mesh.n_middle(obj.solver.mesh.col_to_cell), [], 1);
            dS = X.dX*ones(size(ev_cidx));
            dN = X.N - n_center(ev_cidx); % delta_n
            dZ = zeros(size(ev_cidx)); % delta_z - not used, dSig takes precedence
            dSig = X.Sig - obj.solver.mesh.sig_center(ev_cidx); % delta_sig

            dT = datetime(X.T, 'ConvertFrom', 'datenum');
            M0 = obj.solver.model.get_model(dT, dS, dN, dZ,dSig);

            pars_cell = obj.pars(ev_cidx, :, sol_idx);

            %The following could be coded more elegantly but does not
            %affect computation times significantly (probably
            %matmult/pagemtimes/permute)
            pidx = [0, cumsum(obj.solver.model.npars)];
            
            ip = repmat(1:nx, [np, 1]); ip = ip(:);
            jp = 1:nx*np;
            res = nan([nx, obj.solver.model.ncomponents]);
            v = cell([3,1]); P = cell([3,1]); Mp = cell([3,1]); 
            for dim = 1:obj.solver.model.ncomponents
                v{dim} = M0(:,:,dim)'; v{dim} = v{dim}(:);
                P{dim} = pars_cell(:,(pidx(dim)+1):pidx(dim+1))'; P{dim} = P{dim}(:);
                Mp{dim} = sparse(ip,jp,v{dim});
                res(:,dim) = Mp{dim}*P{dim}; % Compute u = Mx*px;
            end
        end

        function [H, Wl, Zb] = get_H(obj)
            [xvec, yvec] = obj.solver.mesh.xs.sn2xy(zeros(size(obj.decomp.X.y)), obj.decomp.X.y); % convert to x,y
            zvec = obj.bathy.get_bed_elev(xvec, yvec);  % get (time-indep) bathy

            wl = obj.solver.adcp.water_level_object.get_water_level_model(datetime(obj.decomp.X.t, 'ConvertFrom', 'datenum'));


            [Wl, Zb] = ndgrid(wl, zvec);

            H = Wl - Zb;
        end

        function plot_residuals(obj, var_idx)
            nreg = size(obj.p, 2);
            nn = numel(var_idx);

            for fig_idx = 1:3
                m = makefigure(20,3*numel(var_idx));
                if nreg>1 % Compare different vectors
                    t = tiledlayout(nn, nreg, TileSpacing = "tight", Padding = "tight", TileIndexing = "columnmajor");
                else
                    t = tiledlayout('flow', TileSpacing="tight", Padding="tight");
                end

                t.XLabel.String = 'y [m]';
                t.YLabel.String = 'z [m]';
                t.XLabel.Interpreter = 'latex';
                t.YLabel.Interpreter = 'latex';
                var_tit = {'M\vec{p} - \vec{b}', 'C_1\vec{p}',...
                    'C_2\vec{p}', 'C_3\vec{p}', 'C_4\vec{p}', 'C_5\vec{p} - \vec{\zeta}'};
                meas_name = {'m', 'var', 'mse'};
                meas_tit = {'n_s^{-1}\Sigma', 'Var', 'MSE'};
                for col = 1:nreg
                    for row = 1:nn
                        nexttile;
                        obj.plot_residual(col, meas_name{fig_idx}, var_idx{row})

                        lam = {'\mathbf{\lambda}_0', '\mathbf{\lambda}_1', '\mathbf{\lambda}_2'};
                        if row == 1
                            title(['$', meas_tit{fig_idx}, '(',  var_tit{var_idx{row}+1}, ')', ', \hspace{.1cm}  \lambda = ', lam{col}, '$'], 'interpreter', 'latex', 'FontSize', 12);
                        else
                            title(['$', meas_tit{fig_idx}, '(',  var_tit{var_idx{row}+1}, ')', '$'], 'interpreter', 'latex', 'FontSize', 12);
                        end %['$', titles{row}, ', \hspace{.1cm}  \lambda = ', lam{col}, '$']
                        %tit = strcat();
                        %title(tit, 'interpreter', 'latex', 'FontSize', 12);
                        c = colorbar();
                        set(c,'TickLabelInterpreter','latex')
                        colormap(gca, flipud(obj.vel_cmap))
                        c.FontSize = 12;

                        axis tight
                        set(gca, 'XDir','reverse') % Very important

                        set(gca, 'XTick',[])
                        set(gca, 'YTick',[])

                    end
                end
            end

        end

        function plot_solution(obj, plot_names, opts)
            % Wrapper for mesh.plot()
            % Does not support quiver plots
            arguments
                obj
                plot_names (1,:) cell = obj.solver.model.all_names
                opts.sol_idx (1,:) double  = 1:size(obj.p,2)
                opts.representation = "ab"
            end

            if ~isscalar(obj)
                obj.run_method('plot_solution');
                return
            end

            nn = length(plot_names);
            np = sum(obj.solver.model.npars); % Number of parameters in each cell
            m = makefigure(9, 3*nn);
            
            if numel(opts.sol_idx)>1 % Compare different vectors
                t = tiledlayout(m, nn, numel(opts.sol_idx), TileSpacing = "tight", Padding = "tight", TileIndexing = "columnmajor");
            else
                t = tiledlayout(m, 'flow', TileSpacing="tight", Padding="tight");
            end

            t.XLabel.String = 'y';
            t.YLabel.String = 'z';
            t.XLabel.Interpreter = 'latex';
            t.YLabel.Interpreter = 'latex';

            par_idx = obj.get_par_idx(plot_names);
            if strcmp(opts.representation, "Aphi")
                %pa = obj.p2pars();
                [tid_pars, tid_names] = ab2Aphi(obj.pars, obj.solver.model.all_names);
                P = obj.pars2p(tid_pars);
                P = P(:, opts.sol_idx);
                plot_names = tid_names(par_idx);
            else
                P = obj.p(:, opts.sol_idx);
            end

            titles = obj.modify_names(plot_names);
            
            for col = 1:numel(opts.sol_idx)
                for row = 1:nn
                    ax = nexttile(t);
                    amax = max(abs(P(par_idx(row):np:end, numel(opts.sol_idx))), [], 'omitnan') + 1e-5; % For paper, we assume nreg = 3
                    var = P(par_idx(row):np:end, col);
                    hold on
                    [hbed, hwater, hmesh] = obj.solver.mesh.plot(ax,'var', var, 'FixAspectRatio', false);
                    hmesh.LineStyle = 'none';
                    
                    % Colormap
                    % Default: velmap
                    colormap(ax, helpers.cmaps("velmap"))
                    clim([-amax, amax])
                    if contains(titles{row}, 'A')    % linear variable
                        clim([-amax, amax])
                        colormap(ax, helpers.cmaps("Amap"))
                    elseif contains(titles{row}, 'phi')
                        clim([-pi, pi])              % cyclic variable
                        colormap(ax, helpers.cmaps("phimap"))
                    end
                    
                    % Titles
                    lam = {'\mathbf{\lambda}_N', '\mathbf{\lambda}_L', '\mathbf{\lambda}_H'};
                    if 0
                        title(['$', titles{row}, ', \hspace{.1cm}  \lambda = ', lam{col}, '$'], 'interpreter', 'latex', 'FontSize', 12);
                    else
                        title(['$', titles{row}, '$'], 'interpreter', 'latex', 'FontSize', 12);
                    end
                    
                    % Colorbar
                    if col == numel(opts.sol_idx)
                        c = colorbar;
                        c.TickLabelInterpreter = 'latex';
                        c.FontSize = 10;
                        if false
                            if ~contains(titles{row}, '\partial') % Velocities
                                unit = '[m/s]';
                            elseif contains(titles{row}, '\sigma') % Velocities
                                unit = '[m/s]';
                            else % derivatives of velocities in x,y,z directions: m/s/m = 1/s
                                unit = '[1/s]';
                            end



                            pos = get(c,'Position');
                            if row == 1
                                pos1 = pos;
                            end
                            % disp(pos)
                            c.Label.String = unit;
                            c.Label.Interpreter = 'latex';
                            %c.Label.Position(1) = pos1(1) + 3; % to change its position
                            %c.Label.Position(2) = c.Label.Position(2) + .3; % to change its position
                            c.Label.VerticalAlignment = 'middle'; % to change its position
                            c.Label.HorizontalAlignment = 'right';
                            
                            c.Label.Rotation = 270; % to rotate the text
                            
                        end
                    end
                    axis tight
                    set(gca, 'XDir','reverse') % Very important

                    set(gca, 'XTick',[])
                    set(gca, 'YTick',[])

                end
            end
            % Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
            % %[Left Bottom Right Top] spacing
            % NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
            % set(gca, 'Position', NewPos);
        end

        function plot_solution_quiv(obj, names, opts)
            % Wrapper for mesh.plot()
            arguments
                obj
                names (1,:) cell = obj.solver.model.all_names
                opts.sol_idx (1,:) double  = 1:size(obj.p,2)
                opts.representation = "ab"
            end

            if ~isscalar(obj)
                obj.run_method('plot_solution');
                return
            end





%             inp = inputParser;
%             %inp.addOptional('var',[]);
%             inp.addParameter('v', 0, @(x) isscalar(x) && isfinite(x));
%             inp.addParameter('w', 0, @(x) isscalar(x) && isfinite(x));
%             expectedTransform = {'linear','symlog'};
%             inp.addParameter('ArrowScaling',[.1, .1]);
%             inp.addParameter('ArrowTransform','linear', @(x) any(validatestring(x,expectedTransform)));
%             inp.addParameter('ArrowParam', [.9, .9])
%             inp.parse(varargin{:})
%             v = inp.Results.v;
%             w = inp.Results.w;
%             ArrowTransform = inp.Results.ArrowTransform;
%             ArrowScaling = inp.Results.ArrowScaling;
%             ArrowParam = inp.Results.ArrowParam;
            P = obj.p(:, opts.sol_idx);
            nreg = size(P, 2);
            nn = length(names);
            nc = obj.solver.mesh.ncells;
            np = sum(obj.solver.model.npars); % Number of parameters in each cell
            Np = size(P,1); %= nc*np;
            m = makefigure(20, 3*nn);
            if nreg>1 % Compare different vectors
                t = tiledlayout(nn, nreg, TileSpacing = "tight", Padding = "tight", TileIndexing = "columnmajor");
            else
                t = tiledlayout('flow', TileSpacing="tight", Padding="tight");
            end

            t.XLabel.String = 'y [m]';
            t.YLabel.String = 'z [m]';
            t.XLabel.Interpreter = 'latex';
            t.YLabel.Interpreter = 'latex';

            par_idx = obj.get_par_idx(names);
            titles = obj.modify_names(names);
            if strcmp(opts.representation, "Aphi")
                pause
            end
            for col = 1:nreg
                for row = 1:nn
                    ax = nexttile();
                    amax = max(abs(P(par_idx(row):np:Np, 2)), [], 'omitnan') + 1e-5; % For paper, we assume nreg = 3
                    if ~v && ~w
                        var = [P(par_idx(row):np:Np, col), zeros(nc,2)];
                        armax = [0, 0];
                    elseif v && ~w
                        var = [P(par_idx(row):np:Np, col), P(par_idx(row)+np/3:np:Np, col), zeros(nc,1)];
                        armax = [max(abs(P(par_idx(row)+np/3:np:Np, 2)), [], 'omitnan') + 1e-5, 0];
                    elseif w && ~v
                        var = [P(par_idx(row):np:Np, col), zeros(nc,1), P(par_idx(row)+2*np/3:np:Np, col)];
                        armax = [0, max(abs(P(par_idx(row)+2*np/3:np:Np, 2)), [], 'omitnan') + 1e-5];
                    elseif v && w
                        var = [P(par_idx(row):np:Np, col), P(par_idx(row)+np/3:np:Np, col), P(par_idx(row)+2*np/3:np:Np, col)];
                        armax = [max(abs(P(par_idx(row)+np/3:np:Np, 2)), [], 'omitnan') + 1e-5, max(abs(P(par_idx(row)+2*np/3:np:Np, 2)), [], 'omitnan') + 1e-5];
                    end
                    var = obj.arrow_scale(var, ArrowScaling.*armax, ArrowTransform, ArrowParam);
                    %plot(var(:,2:3))
                    %ylim([-armax(1), armax(1)])
                    hold on
                    obj.solver.mesh.plot(ax,'var', var, 'FixAspectRatio', false)
                    if col == nreg
                        if ~contains(titles{row}, 'phi')
                            %amax = max(abs(var(:,1)), [], 'omitnan') + 1e-5;
                            c = colorbar;
                            set(c,'TickLabelInterpreter','latex')
                        else
                            c = colorbar;
                            ylabel(c, 'deg','Rotation',270, 'interpreter', 'latex');

                        end
                    end
                    lam = {'\mathbf{\lambda}_0', '\mathbf{\lambda}_1', '\mathbf{\lambda}_2'};

                    if ~contains(titles{row}, '\partial') % Velocities
                        unit = '[m/s]';
                    elseif contains(titles{row}, '\sigma') % Velocities
                        unit = '[m/s]';
                    else % derivatives of velocities in x,y,z directions: m/s/m = 1/s
                        unit = '[1/s]';
                    end
                    if row == 1
                        title(['$', titles{row}, ', \hspace{.1cm}  \lambda = ', lam{col}, '$'], 'interpreter', 'latex', 'FontSize', 12);
                    else
                        title(['$', titles{row}, '$'], 'interpreter', 'latex', 'FontSize', 12);
                    end

                    if col == nreg
                        pos = get(c,'Position');
                        if row == 1
                            pos1 = pos;
                        end
                        % disp(pos)
                        c.Label.String = unit;
                        c.Label.Interpreter = 'latex';
                        %c.Label.Position(1) = pos1(1) + 3; % to change its position
                        %c.Label.Position(2) = c.Label.Position(2) + .2; % to change its position
                        c.Label.HorizontalAlignment = 'center'; % to change its position
                        c.TickLabelInterpreter = 'latex';
                        c.Label.Rotation = 270; % to rotate the text
                        c.FontSize = 12;
                    end

                    % colormap
                    if ~contains(titles{row}, 'phi')    % linear variable
                        caxis([-amax, amax])
                        %                     colormap(gca, obj.vel_cmap)
                        %                         if col == nreg
                        %                         if amax < 1e-2 % change colorbar ticklabels and ylabel to remove scientific notation
                        %                             c.Ticks = 100*c.Ticks;
                        %                             c.Label.String = [c.Label.String, "$\times 10^{-2}$"];
                        %                         end
                        %                         end
                    else
                        caxis([-180, 180])              % cyclic variable
                        temp = get(gca, 'Children');
                        temp(2).CData =  temp(2).CData*180/pi;
                        %                     colormap(gca, obj.phi_cmap)
                    end

                    axis tight
                    set(gca, 'XDir','reverse') % Very important

                    set(gca, 'XTick',[])
                    set(gca, 'YTick',[])

                end
            end
            % Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
            % %[Left Bottom Right Top] spacing
            % NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
            % set(gca, 'Position', NewPos);
        end

        function ax = plot_mesh(obj, ax, opts)
            arguments
                obj
                ax
                opts.vertical = "z"
            end
            if strcmp(opts.vertical,'z') %dirty
                co = 'z';
                obj.solver.mesh.plot(ax);
            elseif strcmp(opts.vertical,'sig')
                co = '$\sigma$';
                obj.solver.mesh.plot(ax, Sigma = true);
            end 
            xlabel('y')
            ylabel(co, 'Interpreter', 'latex')
            axis tight
            set(gca(), 'XDir','reverse');
            ax.TickLabelInterpreter = 'latex';
            set(gca,'TickLength',[0 0]);
            set(gca, 'XTickLabel', [])
            set(gca, 'YTickLabel', [])
        end

        function CV = cross_validate_single(obj, reg_pars_cell)

            p_train = cell(size(reg_pars_cell));

            if strcmp(obj.opts.cv_mode, 'random')
                nepochs = obj.opts.cv_iter;
            else
                nepochs = 1;
            end
            niter = numel(reg_pars_cell)*nepochs;
            E = cell([numel(reg_pars_cell), 2]);
            CV = cell([numel(reg_pars_cell), 2]);
            i = 0;
            fprintf('Cross-validation percentage: %2.2f percent \n', 100*i/niter)
            for ep = 1:nepochs % randomized iterations: Only meaningful if cv_mode == 'random'
                train_idx = logical(obj.split_dataset());
                test_idx = ~train_idx;

                % Construct training matrix and data
                M0 = obj.M(train_idx, :);
                b0 = obj.b(train_idx);

                M1 = obj.M(test_idx,:);
                b1 = obj.b(test_idx);

                Mp = M0'*M0;

                for rp = 1:numel(reg_pars_cell)
                    regp = reg_pars_cell{rp};
                    p_train{rp}(:, ep) = obj.assemble_solve_single(M0, b0, Mp, regp);
                    i = i+1;
                    fprintf('Cross-validation percentage: %2.2f percent \n', 100*i/niter)
                    E{rp,1}(1, ep) = mean((M1*p_train{rp}(:, ep) - b1).^2); % Generalization error
                    E{rp,2}(1, ep) = mean((M0*p_train{rp}(:, ep) - b0).^2); % Training error
                end
            end
            for rp = 1:numel(reg_pars_cell)
                %                 %p_avg{rp} = mean(p_train{rp}, 2); % k-fold cross-validated estimate
                CV{rp, 1} = mean(E{rp,1}); % Ensemble average
                CV{rp, 2} = mean(E{rp,2}); % Ensemble average
            end
        end


    end

    methods(Access=protected)
        function b_pred = get_b_pred(obj)
            if ~isscalar(obj)
                b_pred = obj.run_method('get_b_pred');
                return
            end
            b_pred = obj.M*obj.p;
        end

        function [res, res_std, res_mean] = get_res(obj)
            if ~isscalar(obj)
                [res, res_std, res_mean] = obj.run_method('get_res');
                return
            end
            res = obj.b - obj.get_b_pred;
            res_std = std(res);
            res_mean = mean(res);
        end

        function [me, vare, mse, e] = get_residual(obj, A, x, b)
            % Input: matrix A, solution x, data b
            % Ax = b + e
            % OLS estimator x will produce mimimal residual norm
            % Biased estimators will produce larger residual norms
            % Output:
            % me = mean(e)
            % vare = sample variance of e, computed as 1/(n-1)*sum(e_i - mean(e))^2
            % mse = mean squared error defined as e^Te/(N-1)
            % e = residual vector

            e = A*x-b;
            [me, vare, mse] = obj.get_mean_var_mse(e);
        end

        function [me, vare, mse] = get_mean_var_mse(obj, e)
            me = mean(e, 'omitnan');
            mse = e'*e./(numel(e)-1);
            if isnan(me)
                mse = nan;
            end
            vare = mse - numel(e)*me^2/(numel(e)-1);
        end

        function idx = eq2mesh(obj, neq)
            % Mapping between equation index and
            % mesh cell index

            % neq(cell_idx) = number of equations within cell_idx
            % idx{cell_idx} = equation indices within cell_idx
            idx = cell([numel(neq),1]);
            cur_idx = 1;
            for cidx = 1:numel(neq)
                idx{cidx} = [cur_idx:(cur_idx + neq(cidx) - 1)]';
                cur_idx = cur_idx + neq(cidx);
            end
        end

        function [me_mesh, vare_mesh, mse_mesh, e_mesh] = get_residual_mesh(obj, A, x, b, neq)
            % Function that obtains the residuals for each mesh cell.
            % Provide the large matrix, solution vector, and rhs
            % Also provide the number of eqs per mesh cell

            % Output: same as get_residual, but now splitted out between
            % mesh cells (added dimension of size mesh.ncells)

            % Implementation using for loop, not very efficient

            me_mesh = nan([obj.solver.mesh.ncells,1]);
            vare_mesh = nan([obj.solver.mesh.ncells,1]);
            mse_mesh = nan([obj.solver.mesh.ncells,1]);
            e_mesh = cell([obj.solver.mesh.ncells,1]);

            idx = obj.eq2mesh(neq);

            [~, ~, ~, e] = obj.get_residual(A, x, b); % Bulk statistics

            for cidx = 1:obj.solver.mesh.ncells
                e_mesh{cidx,1} = e(idx{cidx}, 1);
                [me_mesh(cidx,1), vare_mesh(cidx, 1),  mse_mesh(cidx, 1)] = obj.get_mean_var_mse(e_mesh{cidx,1});
            end
        end

        function gof=get_residuals(obj)
            % Data residuals per mesh:
            for n = 1:size(obj.p,2)
                % Data
                [obj.GOF(n).m, obj.GOF(n).var, obj.GOF(n).mse, obj.GOF(n).e] = obj.get_residual(obj.M, obj.p(:,n), obj.b);
                [obj.GOF(n).mm, obj.GOF(n).varm, obj.GOF(n).msem, obj.GOF(n).em] = obj.get_residual_mesh(obj.M, obj.p(:,n), obj.b, obj.ns);
                % Constraints
                for nc = 1:5
                    [obj.GOF(n).cm{nc}, obj.GOF(n).cvar{nc}, obj.GOF(n).cmse{nc}, obj.GOF(n).ce{nc}]...
                        = obj.get_residual(obj.regularization(nc).C, obj.p(:,n), obj.regularization(nc).rhs);
                    %[obj.GOF(n).cmm{nc}, obj.GOF(n).cvarm{nc}, obj.GOF(n).cmsem{nc}, obj.GOF(n).cem{nc}]...
                    %     = obj.get_residual_mesh(obj.regularization(nc).C, obj.p(:,n), obj.regularization(nc).rhs, obj.regularization.neq{nc});
                end
            end
            gof = obj.GOF;
        end


        function plot_residual(obj, p_idx, meas_name, var_idx)
            % Plots mesh-based residuals with respect do the data and
            % regularization constraints
            % measure_name = subarray of {"m", "var", "mse"}
            % var_idx = 0,1,2,3,4,5 (0: data, 1 - 5: constraints)

            if var_idx == 0
                fi = strcat(meas_name, "m");
                var = obj.GOF(p_idx).(fi);
                %var(var==0) = nan;
                obj.solver.mesh.plot('var', var, 'FixAspectRatio', false)
                amax = max(abs(var), [], 'omitnan');
                %if strmp(meas_name, m)
                clim([-amax, amax])
                %end
            else
                fi = strcat("c", meas_name, "m");
                var = obj.GOF(p_idx).(fi){var_idx};
                amax = max(abs(var), [], 'omitnan');
                %var(var==0) = nan;
                obj.solver.mesh.plot('var', var, 'FixAspectRatio', false)
                clim([-amax, amax])
            end
        end

        function scaled_var = arrow_scale(obj, var, as, at, ap)%, ArrowScaling, ArrowTransform, ArrowParam);
            % Input: N x 3 array, acts on the last two columns
            scaled_var = var;
            if strcmp(at, 'linear')
                scaled_var(:,2:3) = [var(:,2).*as(1), var(:,3).*as(2)];
            elseif strcmp(at, 'symlog')
                scaled_var(:,2:3) = [helpers.symlog(var(:,2), ap(1)).*as(1), helpers.symlog(var(:,3), ap(2)).*as(2)];
            end
        end

        function par_idx = get_par_idx(obj, names_selection)
            par_idx = nan([1,numel(names_selection)]);
            for name_idx = 1:numel(names_selection)
                par_idx(name_idx) = find(strcmp(obj.solver.model.all_names, names_selection{name_idx}));
            end
        end

        function mod_names = modify_names(obj, names_selection)
            mod_names = names_selection;
            for idx = 1:numel(names_selection)
                mod_names{idx} = strrep(mod_names{idx}, 'sig', '\sigma');
                mod_names{idx} = strrep(mod_names{idx}, '^1', '');
                mod_names{idx} = strrep(mod_names{idx}, 'u0', 'u_0');
                mod_names{idx} = strrep(mod_names{idx}, 'v0', 'v_0');
                mod_names{idx} = strrep(mod_names{idx}, 'w0', 'w_0');

                mod_names{idx} = strrep(mod_names{idx}, 'ds', '\partial x');
                mod_names{idx} = strrep(mod_names{idx}, 'dn', '\partial y');
                mod_names{idx} = strrep(mod_names{idx}, 'd\sigma', '\partial \sigma');

                mod_names{idx} = strrep(mod_names{idx}, 'M0', 'M_0');
                mod_names{idx} = strrep(mod_names{idx}, 'M2', 'M_2');
                mod_names{idx} = strrep(mod_names{idx}, 'M4', 'M_4');

                mod_names{idx} = strrep(mod_names{idx}, 'A', '-A');
                mod_names{idx} = strrep(mod_names{idx}, '\phi', '-\phi');
                %mod_names{idx} = strrep(mod_names{idx}, 'M4', 'M_4');
            end
        end

        function training_idx = split_dataset(obj)
            %ci = vertcat(obj.cell_idx{:});
            ci =obj.cell_idx;
            training_idx = ones(size(ci));

            if strcmp(obj.opts.cv_mode, 'none')
                training_idx = ones(size(ci));
            elseif strcmp(obj.opts.cv_mode, 'random')
                tp = obj.opts.training_perc;
                rand0 = rand(size(training_idx));
                training_idx = (rand0 <= tp);
            elseif strcmp(obj.opts.cv_mode, 'omit_cells')
                oc = obj.opts.omit_cells;
                for occ = 1:length(oc)
                    training_idx(ci == oc(occ)) = 0;
                end
            elseif strcmp(obj.opts.cv_mode, 'omit_time') % to be implemented
            end
        end



        function [A, rhs, L] = assemble_single(obj, M, b, Mp, regp)
            A = Mp + regp(1)*obj.solver.regularization(1).Cg+ regp(2)*obj.solver.regularization(2).Cg +...
                regp(3)*obj.solver.regularization(3).Cg + regp(4)*obj.solver.regularization(4).Cg + regp(5)*obj.solver.regularization(5).Cg;
            obj.opts.preconditioner_opts.diagcomp = max(sum(abs(A),2)./diag(A))-2;
            L = ichol(A, obj.opts.preconditioner_opts);
            rhs = M'*b + regp(5)*obj.solver.regularization(5).C'*obj.solver.regularization(5).rhs;
        end

        function p = assemble_solve_single(obj, M, b, Mp, regp)
            % Solve system of eqs
            [A, rhs, L] = assemble_single(obj, M, b, Mp, regp);
            p = obj.solve_single(A, rhs, L); % Matrix of solutions (columns) belonging to regularization parameters regP (rows)
        end

        function p = solve_single(obj, A, rhs, L)
            % Solve system of eqs
            [p, ~, ~, ~] = pcg(A, rhs, obj.opts.pcg_tol, obj.opts.pcg_iter, L, L'); % Matrix of solutions (columns) belonging to regularization parameters regP (rows)
        end

        function rhs = b2rhs(obj, M, b, regp)
            rhs =  M'*b + regp(5)*obj.regularization.C{5}'*obj.regularization.rhs;
        end

        function CV = cross_validation(obj)
            % First loop: cross-validation ensembles
            CV = obj.cross_validate_single(obj.opts.reg_pars, obj.p);
        end

        %function solve_single()

        function CV = cross_validation_analysis(obj)
            reg_pars_cell = reg_pars2cell(obj);
            CV = obj.cross_validate_single(reg_pars_cell, zeros(size(obj.p,1), numel(reg_pars_cell)));
        end

        function SA = local_sensitivity(obj)
            % Perform sensitivity analysis using iid noise ~N(0, sigma^2),
            % only on the state vectors that are estimated in the
            % get_parameters call of Solver.
        end


        function [Pe, P] = local_sensitivity_analysis(obj)
            %             M = dat.M; C1 = dat.C1; C2 = dat.C2; C3 = dat.C3; C4 = dat.C4; C5 = dat.C5; bc = dat.bc; b = dat.b;
            %             Np = size(M,2); % number of parameters
            %             nc = max(dat.cell_idx); % number of cells
            %             np = Np/nc; %number of parameters per cell
            %             opts = dat.opts;
            %             nepochs = opts.cv_iter; %nepochs now serves the role
            %             Mp = M'*M; C1p = C1'*C1; C2p = C2'*C2; C3p = C3'*C3; C4p = C4'*C4; C5p = C5'*C5;
            %             regP = combine_regpars(opts);
            %             pcg_opts = struct('michol','on','type','ict','droptol',1e-3);
            reg_pars_cell = reg_pars2cell(obj);
            if strcmp(obj.opts.generate, 'nullspace') % deprecated
                %     D0 = C3;
                %     for row = 1:size(D0,1) % rescale smoothing matrix back - is this necessary?
                %         D0(row,:) = D0(row,:)/C3(row,row);
                %     end
                %     tic;
                %                 NS = generate_null_intersection({C1, C2, C3, C4, C5});
                %     to = toc;
                %                 fprintf('Finished calculating intersection of null spaces after %2.2f s \n', to)

                %                 p = NS*randn(size(NS,2), nepochs);
            elseif strcmp(obj.opts.generate, 'local')
                Mp = obj.M'*obj.M;
                p0 = nan([size(obj.p,1),numel(reg_pars_cell)]);
                for rp = 1:numel(reg_pars_cell) % generate solution for all regularization parameters, using all data
                    regp = reg_pars_cell{rp};
                    p0(:, rp) = obj.assemble_solve_single(obj.M, obj.b, Mp, regp, zeros([size(obj.p,1),1]));
                    %p0: Solution state vectors belonging to all
                    %regularization parameters in the sensitivity analysis.
                    fprintf('Ensemble generation percentage: %2.2f percent \n', 100*rp/numel(reg_pars_cell))
                end
                %p0 = mean(p,2); % Selection of weighted parameter vectors across many reg.pars.
            end


            W0 = rand([numel(reg_pars_cell), obj.opts.ens_size]);
            W = W0/diag(sum(W0,1)); %Scale all columns
            P = p0*W; % True vectors: columns of P

            stdn = obj.opts.get_noise();

            B = obj.M*P; % Unperturbed data: ns x nP
            niter = numel(reg_pars_cell)*numel(stdn)*obj.opts.ens_size*obj.opts.sa_iter; % total number of pcg runs
            fprintf('Total number of iterations will be: %i \n', niter)
            i=0;
            Mp = obj.M'*obj.M;
            Pe = cell([numel(reg_pars_cell), numel(stdn)]);
            %E = cell([numel(reg_pars_cell), numel(stdn)]);
            for rp = 1:numel(reg_pars_cell)
                regp = reg_pars_cell{rp};
                [A, ~, L] = obj.assemble_single(obj.M, obj.b, Mp, regp);
                for nn = 1:numel(stdn)
                    for it = 1:obj.opts.sa_iter % Bootstrapping random iterations
                        Bp = B + stdn(nn)*randn(size(B));
                        for ep = 1:obj.opts.ens_size
                            % TODO apply variation in sa_iter: Noise term must change.
                            i = i+1;
                            % Perturbed measurements -> perhaps move to inner loop.
                            rhs = obj.b2rhs(obj.M, Bp(:,ep), regp); % Estimate different base vector per iteration.
                            fprintf('Sensitivity analysis: %2.2f percent \n', 100*i/niter)
                            Pe{rp, nn}(:, ep, it) = obj.solve_single(A, rhs, L, P(:,ep));
                        end
                    end
                end
            end
        end

        function E = extract_sensitivity_data(obj, Pe, P)
            for rp = 1:size(Pe, 1)
                for nn = 1:size(Pe, 2)
                    for it = 1:obj.opts.sa_iter % Bootstrapping random iterations
                        for ep = 1:obj.opts.ens_size
                            E{rp, nn, 1}(ep, it) = mean((P(:, ep) - Pe{rp, nn}(:, ep, it)).^2); %MSE
                            E{rp, nn, 2}(ep, it) = mean((P(:, ep) - Pe{rp, nn}(:, ep, it)).^2);
                        end
                    end
                end
            end
        end


        function reg_pars_cell = reg_pars2cell(obj)
            reg_pars_sens_vec = obj.opts.vectorize_reg_pars();
            reg_pars_cell = cell(size(reg_pars_sens_vec,1), 1);
            for rp = 1:size(reg_pars_sens_vec,1)
                reg_pars_cell{rp} = reg_pars_sens_vec(rp,:);
            end
        end

        function fig = plot_contourf_fig(obj, var, xvar, yvar, scaled, tit)
            fig = makefigure(24,12);
            t = tiledlayout('flow', TileSpacing="tight", Padding="tight");
            %             t = tiledlayout(2,3, TileSpacing="tight", Padding="tight");
            %s=nexttile;

            t.XLabel.String = '$\lambda_c$';
            t.YLabel.String = '$\lambda_g$';
            t.XLabel.Interpreter = 'latex';
            t.YLabel.Interpreter = 'latex';
            for idx = 1:numel(var)

                ax{idx} = nexttile;
                plot_contourf_ax(obj, ax{idx}, var{idx}, xvar, yvar, scaled{idx}, tit{idx});
            end
        end

        function ax = plot_contourf_ax(obj, ax, var, xvar, yvar, scaled, tit)

            if scaled
                var = helpers.matdivrobust(var, var(1,1));
                caxis([1-max(max(abs(var-1))), 1+max(max(abs(var-1)))]);
            end
            hold on
            contourf(xvar, yvar, var, 50, 'LineStyle','none')
            contour(xvar, yvar, var, [1,1], 'LineColor','k')
            hold off
            %xlabel('lc'); ylabel('lg')
            colormap(gca, flipud(brewermap(50, 'RdBu')));

            title(tit, 'interpreter', 'latex')
            idxt = unique([1, round(size(xvar,2)/4), round(size(xvar,2)/2), round(3*size(xvar,2)/4), round(size(xvar,2))]);
            idyt = unique([1, round(size(yvar,1)/4), round(size(yvar,1)/2), round(3*size(yvar,1)/4), round(size(yvar,1))]);
            ax.XTick = xvar(1,idxt);
            ax.YTick = yvar(idyt,1);
            ax.XTickLabel = round(helpers.symexp(xvar(1,idxt), obj.opts.res_near_zero(1)), 1, 'significant');
            ax.YTickLabel = round(helpers.symexp(yvar(idyt,1), obj.opts.res_near_zero(3)), 1, 'significant');
            ax.TickLabelInterpreter = 'latex';
            c = colorbar();
            c.TickLabelInterpreter = 'latex';
        end


        function plot_cross_validation_analysis(obj, CV, scaled)
            [CV, lc, lg] = prep_sens_plot(obj, CV);

            fig = plot_contourf_fig(obj, CV, lc, lg, {1,1}, {'Generalization error', 'Training error'});


        end

        function [plot_var, lc, lg] = prep_sens_plot(obj, SA)
            plot_var = cell([size(SA,2), 1]);
            for idx = 1:size(SA,2) % number of metrics to be plotted
                plot_var{idx,1} = reshape(cell2mat(SA(:,idx)),...
                    [obj.opts.reg_iter(3), obj.opts.reg_iter(1)]);
                % CV contains metrics belonging to different reg-pars
            end
            if strcmp(obj.opts.reg_vary, 'coupled')
                [lc, lg] = meshgrid(helpers.symlog(obj.opts.reg_pars_sens{1}, obj.opts.res_near_zero(1)),...
                    helpers.symlog(obj.opts.reg_pars_sens{3}, obj.opts.res_near_zero(3)));
            else
                warning("Use SolverOptions.reg_vary = 'coupled' to perform sensitivity analyses.")
            end
            %             plot_var
            %             lc
            %             lg
        end

        function pars = p2pars(obj,p)
            [np, nsols] = size(p);
            ncells = obj.solver.mesh.ncells;
            npars = np/ncells;
            pars = reshape(obj.p, npars, ncells, nsols);
            pars = permute(pars, [2 1 3]);
        end

        function p = pars2p(obj, pars)
            p = nan(size(obj.p));
            for i=1:size(pars,3)
                p(:,i) = reshape(pars(:,:,i)', 1, [])';
            end
        end

        function nullSpace = generate_null_intersection(obj, mat_cell)
            % Very slow function that takes cell array of (sparse) matrices and computes
            % intersection of their nullspaces, with basis vectors placed in columns of nullSpace
            cur_null = null(full(mat_cell{1}));
            for m = 2:length(mat_cell)
                next_null = null(full(mat_cell{m}));
                pre_null = null([cur_null, -next_null]);
                cur_null = cur_null*pre_null(1:size(cur_null,2),:);
            end
            nullSpace = cur_null;
        end
    end
end
