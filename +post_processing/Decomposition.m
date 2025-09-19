classdef Decomposition < handle & helpers.ArraySupport
    % Decomposition of functions on time-varying cross-sectional domains
    % Low-level class that decomposes function F: Omega \subset R^3 -> R
    % where Omega is a time-varying domain
    properties
        X struct % Coordinates (t, y, sig, T, Y, Sig)

        H double % Local water depth, depending on (t,y), artificially extended to 3D variable by repmat in the sigma direction for convenience

        wl double

        zb double

        %DFf struct % Orthogonal components of the function to be decomposed
    end
    properties(Dependent)
        A
        B
        sz
    end
    methods
        function obj = Decomposition(varargin)
            obj = obj@helpers.ArraySupport(varargin{:})
        end

        function A = get.A(obj)
            A = obj.wi(obj.H); % Nt by 1 vector
            A = squeeze(A(:,1,1));
        end

        function B = get.B(obj)
            B = abs((obj.X.y(end) - obj.X.y(1)));
        end

        function sz = get.sz(obj)
            sz = size(obj.X.T);
        end

        % Elementary functions - trapz or mean
        function tif = ti(obj, f)
            tif = trapz(obj.X.t, f, 1);
        end

        function wif = wi(obj, f)
            wif = trapz(obj.X.y, f, 2);
        end

        function dif = di(obj, f)
            dif = trapz(obj.X.sig, f, 3);
        end

        % Thickness-weighed functions %TODO look in detail.
        function f = repmat3(obj, f)
            % Reshape matrix to size sz , repeating mean values if needed
            rep = [1 1 1];
            if numel(size(f)) == 2
                rep(3) = obj.sz(3);
            end
            rep([size(f)==1]) = obj.sz([size(f)==1]);
            if ~isempty(rep)
                f = repmat(f, rep);
            end
        end


        function [twaf1, twaf2] = twa(obj, f)
            % Taking thickness-weighed average over t
            % Could be vectorized but for now I prefer readability
            % Returns squeezed AND full array

            [f, sf] = obj.handle_input_dim(f);

            if sf(2) == 1 % Laterally averaged variable
                twaf0 = obj.ti(obj.A.*f)./obj.ti(obj.A);
            else
                twaf0 = obj.ti(obj.H.*f)./obj.ti(obj.H); % Just a factor B difference. Yes.
            end
            twaf1 = twaf0(1, 1:sf(2), 1:sf(3));
            twaf2 = obj.repmat3(twaf0);
        end

        function [laf1, laf2] = la(obj,f)
            % Returns squeezed AND full array
            [f, sf] = obj.handle_input_dim(f);
            if sf(1) == 1 % Time-averaged variable
                laf0 = obj.wi(obj.ti(obj.H).*f)./obj.ti(obj.A);
            else
                laf0 = obj.wi(obj.H.*f)./obj.wi(obj.H);
            end
            laf1 = laf0(1:sf(1), 1, 1:sf(3));
            laf2 = obj.repmat3(laf0);
        end

        function [vaf1, vaf2] = da(obj,f)
            % Returns squeezed AND full array
            [f, sf] = obj.handle_input_dim(f);
            vaf0 = obj.di(f)/(obj.X.sig(end)-obj.X.sig(1));
            vaf1 = vaf0(1:sf(1), 1:sf(2), 1);
            vaf2 = obj.repmat3(vaf0);
        end

        function [f, sf] = handle_input_dim(obj,f)
            sf = size(f); % Get original size of input data
            if size(f,3)==1
                sf = [sf 1];
            end
            f = obj.repmat3(f);
        end

        function plot_components(obj, AF, q)
            figure('units','normalized','outerposition',[0 0 1 1])
            for i=1:numel(AF)
                subplot(2,4,i)
                obj.plot_component(AF{i})
                %                 ncolor=100;
                colormap(gca, helpers.cmaps(q));
            end
        end

        function plot_component(obj, af) %TODO change
            sf = size(af);
            if size(af,3)==1
                sf = [sf 1];
            end
            if sf(1) == 1 % t-independent
                if sf(2) == 1 % y-independent
                    if sf(3) == 1 % f0
                        bar(af)                                      %1
                        xticklabels("f0")
                        title("F0")
                    else % only sig-dependent
                        v = squeeze(af);
                        plot(v, obj.X.sig)    %2
                        xlabel("f(sig)")
                        ylabel("sig")
                        title("F(sig)")
                        xlim([-max(abs(v), [], "all")-1e-6, max(abs(v), [], "all")+1e-6])
                    end
                else % y-dependent
                    if sf(3) == 1 % only y-dependentcma
                        v = squeeze(af);
                        plot(obj.X.y, v)     %3
                        xlabel("y")
                        ylabel("f(y)")
                        title("F(y)")
                        ylim([-max(abs(v), [], "all")-1e-6, max(abs(v), [], "all")+1e-6])
                    else % y and sig dependent
                        v=squeeze(af)';
                        contourf(obj.X.y, obj.X.sig, v)                        %4
                        xlabel("y")
                        ylabel("sig")
                        title("F(y,sig)")
                        colorbar
                        clim([-max(abs(v), [], "all")-1e-6, max(abs(v), [], "all")+1e-6])
                    end
                end
            else % t-dependent
                if sf(2) == 1 % y-independent
                    if sf(3) == 1 % only t-dependent
                        v = squeeze(af)';
                        plot(obj.X.t, v)
                        xlabel("t")
                        ylabel("f(t)")%5
                        title("F(t)")
                        ylim([-max(abs(v), [], "all")-1e-6, max(abs(v), [], "all")+1e-6])
                    else % only (t,sig)-dependent
                        v =  squeeze(af)';
                        contourf(obj.X.t, obj.X.sig, v)
                        xlabel("t")
                        ylabel("sig")                        %6
                        colorbar
                        title("F(t,sig)")
                        clim([-max(abs(v), [], "all")-1e-6, max(abs(v), [], "all")+1e-6])
                    end
                else % y-dependent
                    if sf(3) == 1 % only (t,y)-dependent
                        v = squeeze(af)';
                        contourf(obj.X.t, obj.X.y, v)
                        xlabel("t")
                        ylabel("y")%7
                        title("F(t,y)")
                        colorbar
                        clim([-max(abs(v), [], "all")-1e-6, max(abs(v), [], "all")+1e-6])
                    else % t and y and sig dependent
                        histogram(af(:))
                        title("F(t,y,sig)")
                        % disp('Can only display at most 2D variables')
                    end
                end
            end
            %             colormap(helpers.cmaps(q))
        end


        function val = extract_inst_vlim(obj, af, vlim)

            [~, idx(1)] = min(abs(obj.X.sig - vlim(1)));
            [~, idx(2)] = min(abs(obj.X.sig - vlim(2)));

            val = 1./(obj.X.sig(idx(2)) - obj.X.sig(idx(1))).*trapz(obj.X.sig(idx(1):idx(2)), af(:,:,idx(1):idx(2)), 3);

        end

        function val = extract_inst_llim(obj, af, llim)
            % af: function
            % llim: 1x2 array with elements increasing in the interval
            % [-.5, .5].

            % val: Array Nt x 1 x Nsig

            [~, idx(1)] = min(abs(obj.X.y - llim(1).*obj.B));
            [~, idx(2)] = min(abs(obj.X.y - llim(2).*obj.B));

            afH = trapz(obj.X.y(idx(1):idx(2)), obj.H(:, idx(1):idx(2), :), 2);

            val = 1./afH.*trapz(obj.X.y(idx(1):idx(2)), obj.H(:, idx(1):idx(2), :).*af(:, idx(1):idx(2), :), 2);

        end

        function val = extract_on_tgrid(obj, af, tgrid)
            val = af(tgrid, :, :);
        end

        function [f, ff] = avg_comp(obj, f, avg_op)
            % Arbitrary composition of averaging operators
            sf0 = size(f);
            if size(f,3)==1
                sf0 = [sf0 1];
            end
            %f0 = f;
            ff = obj.repmat3(f);
            if numel(avg_op) > 0
                for idx = 1:numel(avg_op)
                    if isa(avg_op{idx}, 'function_handle')
                        [f, ~] = avg_op{idx}(f);
                        [~, ff] = avg_op{idx}(ff);
                    end
                end
            end
        end

        function [f, ff] = avg_0(obj, f)
            avg_op = {@obj.twa, @obj.la, @obj.da};
            [f,ff] = obj.avg_comp(f, avg_op);
        end

        function [names, ap] = get_names(obj, vname)
            ap = ["_0", "}_t^t", "_y^y", "_z^z", "_{yz}", "_{tz}", "_{ty}", "_{tyz}"];
            op = ["", "\bar{", "\bar{", "\bar{", "\hat{", "\underline{", "[", ""];
            for idx = 1:numel(ap)
                names{idx, :} = ['$',  op(idx), vname, ap(idx), '$'];
            end
        end


        function [DF, AF] = decompose_function(obj, f)
            % DF: Orthogonal basis functions for constructing f
            % Averaged components, non-orthogonal
            % f_0
            % bar(f)(t)
            % bar(f)(y)
            % bar(f)(sig)
            % twa(f)
            % la(f)
            % da(f)
            % full(f)

            AF{8} = f;

            AF{7} = obj.avg_comp(f, {@obj.da});     % da(f) -> depth averaged
            AF{6} = obj.avg_comp(f, {@obj.la});     % la(f) -> laterally averaged
            AF{5} = obj.avg_comp(f, {@obj.twa});    % twa(f) -> tidally averaged

            AF{4} = obj.avg_comp(AF{5}, {@obj.la}); % la(twa(f)) -> residual depth profile
            AF{3} = obj.avg_comp(AF{5}, {@obj.da}); % da(twa(f)) -> residual lateral profile
            AF{2} = obj.avg_comp(AF{7}, {@obj.la}); % la(da(f)) ->  CSA temporal profile

            AF{1} = obj.avg_comp(AF{3}, {@obj.la}); % f_0

            % Decomposition terms: mutually orthogonal w.r.t. ()_0
            % f_0
            % bar(f)_t(t)
            % bar(f)_y(y)
            % bar(f)_sig(sig)
            % twa(f)_(y sig) (y, sig)
            % la(f)_(t sig) (y, sig)
            % da(f)_(t y) (y, sig)
            % full(f)

            DF{1} = AF{1};

            DF{2} = AF{2} - DF{1};
            DF{3} = AF{3} - DF{1};
            DF{4} = AF{4} - DF{1};

            DF{5} = AF{5} - DF{1}-DF{3}-DF{4};
            DF{6} = AF{6} - DF{1}-DF{2}-DF{4};
            DF{7} = AF{7} - DF{1}-DF{2}-DF{3};

            DF{8} = AF{8} - DF{1}-DF{2}-DF{3}-DF{4}-DF{5}-DF{6}-DF{7};
        end

        function prod_idx = prod_idx(obj)
            % input: cell of indices
            % output: cell of indices and their swapped counterparts
            prod_idx_pre = {[0,1], [2,6], [3,5], [4,7];...
                            [0,2], [1,6], [3,4], [5,7];...
                            [0,3], [1,5], [2,4], [6,7];...
                            [0,4], [1,7], [2,3], [5,6];...
                            [0,5], [1,3], [2,7], [4,6];...
                            [0,6], [1,2], [3,7], [4,5];...
                            [0,7], [1,4], [2,5], [3,6]};

            prod_idx_pre = sym_idx(prod_idx_pre);
            prod_idx = vertcat({[0,0], [1,1], [2,2], [3,3], [4,4], [5,5], [6,6], [7,7]}, prod_idx_pre);
        end

        function avg_op_cell = avg_op_cell(obj)
            avg_op_cell = {{@obj.avg_0};
                {@obj.da, @obj.la};
                {@obj.da, @obj.twa};
                {@obj.twa, @obj.la};
                {@obj.twa};
                {@obj.la};
                {@obj.da};
                {}};
        end

        function USF = get_prod_components_all(obj, DF1, DF2)
            % Get all functions AND relevant correlations constituting ALL
            % transport terms

            % Returns 8x8x3 cell, (:,:,1) u and (:,:,2) s            
            % (:,:,3) is the net contribution to the transport / flux (not
            % neccesarily net residual) -> i.e. averaged over the
            % dimension(s) at hand. i.e. second row = bar{fg}_t^t -> 8 terms averaged
            % over y and sigma.

            % (:,:,4) is square root of the norm of the contribution to the transport.

            % (:,:,5) is the name

            %             SF = decompose_product_inst(obj, DF1, DF2);

            avg_op_cell = obj.avg_op_cell();
            prod_idx = obj.prod_idx();

            USF = cell([8,8,5]);

            for fidx = 1:8
                for cidx = 1:8
                    USF{fidx, cidx, 1} = DF1{prod_idx{fidx, cidx}(1)+1};
                    USF{fidx, cidx, 2} = DF2{prod_idx{fidx, cidx}(2)+1};
                    USF{fidx, cidx, 3} = obj.avg_comp(USF{fidx, cidx, 1}.*USF{fidx, cidx, 2}, avg_op_cell{fidx});
                    USF{fidx, cidx, 4} = obj.avg_0(USF{fidx, cidx, 3}.^2); % square of norm of component
                    USF{fidx, cidx, 5} = ['f', num2str(prod_idx{fidx, cidx}(1)), 'g', num2str(prod_idx{fidx, cidx}(2))];
                end
            end % These are all possible terms in our research!
        end

        function SF = decompose_product_inst(obj, DF1, DF2)
            SF = cell(numel(DF1), numel(DF2));
            for i = 1:numel(DF1)
                for j = 1:numel(DF2)
                    SF{i,j}= DF1{i}.*DF2{j}; % Since we have proven commutativity
                end
            end
        end


        function SF = decompose_product_diag(obj, DF1, DF2)
            SF = zeros(numel(DF1), 1);
            for i = 1:numel(DF1)
                SF(i, 1) = obj.avg_0(DF1{i}.*DF2{i}); % Since we have proven commutativity
            end
        end

        %         function SF = decompose_product_full(obj, DF1, DF2)
        %             SF = cellfun(@obj.avg_0, cellfun(@times, repmat(DF1', [1,8]), repmat(DF2, [8,1]), UniformOutput=false));
        %         end

        function SF = decompose_product_full(obj, DF1, DF2)
            SF = zeros(numel(DF1), numel(DF2));
            for i = 1:numel(DF1)
                for j = 1:numel(DF2)
                    SF(i, j) = obj.avg_0(DF1{i}.*DF2{j}); % Since we have proven commutativity
                end
            end
        end

        function SF = get_norms(obj, DF1, DF2, avg_op)
            SF = obj.avg_op_all(DF1, DF2, avg_op);
            %             all_op = {@obj.twa, @obj.la, @obj.da};
            for i = 1:numel(DF1)
                for j = 1:numel(DF2)
                    SF{i, j} = obj.avg_0(SF{i, j}.^2);
                end
            end
        end



        function plotC(obj, C)
            figure('units','normalized','outerposition',[0 0 1 1])
            imagesc(C)
            colorbar
            pos=get(gca,'position');
            [rows,cols]=size(C);
            width= 1.05*(pos(3)-pos(1))/(cols-1);
            height = (pos(4)-pos(2))/(rows-1);
            % %create textbox annotations
            for i=1:cols
                for j=1:rows
                    annotation('textbox',[pos(1)+width*(i-1),pos(2)+height*(rows-j),width,height], ...
                        'string',num2str(C(i,j)),'LineStyle','none','HorizontalAlignment','center',...
                        'VerticalAlignment','middle');
                end
            end
            xticklabels(["$u_0$", "$\bar{u}_t^t$", "$\bar{u}_y^y$", "$\bar{u}_\sigma^\sigma$",...
                "$\widehat{u}_{y\sigma}^{y\sigma}$", "$\underline{u}_{t\sigma}^{t\sigma}$", "$[u]_{ty}^{ty}$", "$u_{t y\sigma}$"]);
            yticklabels(["$s_0$", "$\bar{s}_t^t$", "$\bar{s}_y^y$", "$\bar{s}_\sigma^\sigma$",...
                "$\widehat{s}_{y\sigma}^{y\sigma}$", "$\underline{s}_{t\sigma}^{t\sigma}$", "$[s]_{ty}^{ty}$", "$s_{t y\sigma}$"]);

        end

        function plotUSF(obj, USF, fidx)
            %figure('units','normalized','outerposition',[0 0 1 1])

            %data
            cidx = 1:8;
            dat = cell2mat(squeeze(USF(fidx+1,:,4)));

            %plotting
            fidxp = 1:numel(fidx);
            [Fidx, Cidx] = meshgrid(fidxp, cidx);            
            imagesc(cidx, fidxp, dat)
            colorbar;
            It = Fidx'; Jt = Cidx';
            text(Jt(:), It(:), squeeze(USF(fidx+1,:,5)), 'HorizontalAlignment', 'Center')
            yticks(cidx)
            xticklabels(cidx)
            yticklabels(fidx)
        end
    end
end

%         function [op, idx] = get_bar(obj, dimension)
%             switch dimension
%                 case "t"
%                     op = {@obj.la, @obj.da};
%                     idx = {[0,1], [1,0], [2,6], [6,2], [3,5], [5,3], [4,7], [7,4] };
% %                     idx
%                 case "y"
%                     op = {@obj.twa, @obj.da};
%                     idx = {[0,2], [2,0], [1,6], [6,1], [3,4], [4,3], [5,7], [7,5]};
%                 case "sig"
%                     op = {@obj.twa, @obj.la};
%                     idx = {[0,3],[3,0],  [1,5], [5,1], [2,4],[4,2], [6,7], [7,6]};
%                 otherwise
%                     op = {};
%                     idx = {};
%             end
%         end
%
%         function [op, idx] = get_avg_single(obj, dimension)
%             switch dimension
%                 case "t"
%                     op = @obj.twa;
%                     idx = {[0,1], [1,0], [2,6], [6,2], [3,5], [5,3], [4,7], [7,4] };
% %                     idx
%                 case "y"
%                     op = @obj.la;
%                     idx = {[0,2], [2,0], [1,6], [6,1], [3,4], [4,3], [5,7], [7,5]};
%                 case "sig"
%                     op = @obj.da;
%                     idx = {[0,3],[3,0],  [1,5], [5,1], [2,4],[4,2], [6,7], [7,6]};
%                 otherwise
%                     op = {};
%                     idx = {};
%             end
%         end


%         function [f,ff] = avg_single(obj, f, dimension)
%             avg_op = obj.get_avg_single(dimension);
%             [f, ff] = obj.avg_comp(f, avg_op);
%         end
%
%         function [f,ff] = bar(obj, f, dimension)
%             avg_op = obj.get_bar(dimension);
%             [f, ff] = obj.avg_comp(f, avg_op);
%         end
%
%         function [f,ff] = bar_t(obj, f)
%             avg_op = obj.get_bar("t");
%             [f, ff] = obj.avg_comp(f, avg_op);
%         end
%
%         function [f,ff] = bar_y(obj, f)
%             avg_op = obj.get_bar("y");
%             [f, ff] = obj.avg_comp(f, avg_op);
%         end
%
%         function [f,ff] = bar_sig(obj, f)
%             avg_op = obj.get_bar("sig");
%             [f, ff] = obj.avg_comp(f, avg_op);
%         end
%



%
%         function SF = avg_op_all(obj, DF1, DF2, avg_op_str)
%             SF = decompose_product_inst(obj, DF1, DF2); %multiply funs to obtain a 8x8 matrix
%             if numel(avg_op_str) >0
%                 for aop = 1:numel(avg_op_str)
%                     avg_op{aop} = obj.get_avg_single(avg_op_str{aop});
%                 end
%             else
%                 avg_op = {};
%             end
%             for i = 1:numel(DF1)
%                 for j = 1:numel(DF2)
%                     SF{i, j} = obj.avg_comp(DF1{i}.*DF2{j}, avg_op); % Since we have proven commutativity
%                 end
%             end
%         end
%
%         function [SF0, SFd] = decompose_product_bar(obj, dimension, DF1, DF2)
%             % Yields elements constituting the
%             % - dimension = "t": tidal salt transport, i.e.
%             % la(da(us))(t) = (us)_0 + (us)_t^t
%
%             % Residual CSA contributions:
%
%             [~, idx] = obj.get_bar(dimension); % indexes of terms yielding net transport
%                 %         obj.bar(    )
%             for i = 1:numel(idx) % 8x1 cell
%                 idxm = idx{i} + 1; % matlab indices
% %                 idxms =
%                 SF0{i} = obj.bar(DF1{i}.*DF2{i}, dimension);
%                 SFd{i} = obj.bar(DF1{idxm(1)}.*DF2{idxm(2)}, dimension); % because of symmetry of indices
%             end
%         end
