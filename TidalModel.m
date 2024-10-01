classdef TidalModel < DataModel
    % Tidal model of the data
    %
    %   TidalModel properties:
    %   constituents
    %
    %   TidalModel methods:
    %   get_tidal_pars - compute the tidal amplitude and phase
    %
    %   see also: DataModel, TaylorModel
    properties(Constant)
        supported_constituents = {'M2','S2','N2','K1','M4','O1','M6',...
            'MK3','S4','MN4'}
    end
    properties
        % TidalModel/constituents constituents to fit data with
        %
        %   1 x n_constituents cell array of tidal constituent names.
        %
        %   see also: TidalModel, get_tidal_pars
        constituents = {};
    end

    properties(Dependent, SetAccess = protected)
        % TidalModel/constituents constituents to fit data with
        %
        %   MxN array defining tidal constituents to be included in the model
        %   for data to be fitted. Every row is a different component, e.g. for
        %   velocity: row 1 is x component, row 2 is y component and row 3 is z
        %   component. For a scalar only one row is provided.
        %   Columns represent different constituents. Zero or NaN values are
        %   skipped. For every given value, two model parameters are fitted, which
        %   are the coefficients for the cos and sin functions. From those amplitude
        %   and phases are computed. A residual component is always included.
        %
        %   see also: TidalModel, get_tidal_pars
        periods

        nconstituents
    end



    methods
        function val = get.nconstituents(obj)
            val = obj.get_nconstituents;
        end
        function set.constituents(obj,val)
            assert(ischar(val) || iscellstr(val),...
                'Give constituents as char or cell of chars') %#ok<ISCLSTR>
            if ischar(val)
                val = {val};
            end
            assert(all(ismember(val, obj.supported_constituents)),...
                'Not all given constituents are supported')
            obj.constituents = val;
        end
        function val = find_par(obj, options)
            arguments(Input)
                obj
                options.component(1,:) string ...
                    {mustBeValidComponent(obj, options.component)} =...
                    obj.get_component_names
                options.constituent(1,:) string ...
                    {mustBeValidConstituent(obj, options.constituent)} =...
                    [{'M0'}, obj.supported_constituents]
                options.sincos(1,:) string ...
                    {mustBeSinOrCos(obj, options.sincos)} = ["sin" "cos"]
            end
            all_comps = obj.get_component_names;
            
            comp = options.component;
            const = options.constituent;
            sincos_offset = [0 1];
            sincos = ismember({'cos', 'sin'},options.sincos);
            sincos = sincos_offset(sincos);
            val = false(sum(obj.get_npars_tid),1);
            pos = 1;
            for cc = 1:obj.ncomponents
                comp_name = all_comps(cc);
                if ismember(comp_name, comp) &&...
                    ismember({'M0'}, const)
                    val(pos) = true;
                end
                pos = pos + 1;

                for ct = 1:numel(obj.constituents)
                    const_name = obj.constituents(ct);
                    if ismember(comp_name, comp) &&...
                        ismember(const_name, const)
                        val(pos+sincos) = true;
                    end
                    pos = pos + 2;
                end
            end
        end
        function M = get_model(obj, time, ~, ~, ~, ~)
            % This model fits the following parameters to the velocity
            % within each cell:
            % u = u_0 + sum_n (a_n cos(2pi/T_n * t) + b_n sin(2pi/T_n * t))
            % where n loops over all entered constituents (subtidal
            % constituent is always present)
            % Input:
            % d_time (seconds)
            % Output:

            npars = obj.get_npars_tid;
            ncomp = obj.ncomponents;
            assert(isdatetime(time), 'Enter time vector in datetime format.')
            % As datenum is in days, convert to seconds
%             d_t = convertTo(d_time,"datenum")*24*3600; %seconds
            % Edit 12-8-24 : converted back to datenum.
            secday = 3600*24;
            tn = datenum(time)*secday;
%             tn = seconds(time - time(1));
            max_pars = max(npars);
            M = nan(numel(time), max_pars, ncomp);
            M(:,1,:) = 1; %residual
            for c_comp = 1:ncomp
                n_const = sum(isfinite(obj.periods(c_comp,:)) &...
                    obj.periods(c_comp,:)~=0);
                for c_const = 1:n_const
                    M(:,2*c_const, c_comp) = cos(...
                        2*pi/obj.periods(c_comp,c_const)*tn);
                    M(:,2*c_const + 1, c_comp) = sin(...
                        2*pi/obj.periods(c_comp,c_const)*tn);
                end
            end
        end

        function [pars_h, cov_pars_h] = get_tidal_pars(obj, pars, cov_pars)
            % Compute tidal amplitude and phase from sin and cos coefficients
            %
            %   [pars_h, cov_pars_h] = get_tidal_pars(obj, pars, cov_pars) computes the
            %   amplitude and phases based on the sin and cos coefficients:
            %   u = A cos(omega t - phi) = a cos(omega t) + b sin(omega t) with:
            %   A = sqrt(a^2 + b^2)
            %   phi = arctan(b / a)
            %
            %   pars is an array holding the tidal model parameters as produced with
            %   the Solver class.
            %
            %   cov_pars holds the covariance matrix of the model parameters.
            %
            %   see also: TidalModel, Solver

            npars = obj.get_npars_tid;
            subtidal_idx = [1, npars(1) + 1, sum(npars(1:2)) + 1];
            pars_h = zeros(size(pars));
            jac = zeros(size(cov_pars)); % jacobian of the transformation (subtidal independent of other under this transformations)
            idx = 1;
            while idx <= size(pars_h,2)
                if ismember(idx, subtidal_idx)
                    pars_h(:,idx) = pars(:,idx);
                    jac(:,idx,idx) = 1;
                    idx = idx + 1;
                else
                    pars_h(:,idx) = sqrt(pars(:,idx).^2 + pars(:,idx+1).^2); % Amplitude of constituent
                    pars_h(:,idx+1) = atan(pars(:,idx+1)./pars(:,idx));     % Phase of constituent
                    jac(:,idx,idx) = pars(:,idx)./pars_h(:,idx);            % d(sqrt(a^2 + b^2))/da = a/A
                    jac(:,idx, idx+1) = pars(:,idx+1)./pars_h(:,idx);       % d(sqrt(a^2 + b^2))/db = b/A
                    jac(:,idx+1,idx) = -pars(:,idx+1)./pars_h(:,idx).^2;    % d(atan(b/a))/da = -b/A^2
                    jac(:,idx+1,idx+1) = pars(:,idx)./pars_h(:,idx).^2;     % d(atan(b/a))/db = a/A^2
                    idx = idx + 2;
                end
            end
            cov_pars_h = helpers.matmult(cov_pars,permute(jac,[1,3,2]),2,3);   % cov_pars * J'
            cov_pars_h = helpers.matmult(jac,cov_pars_h,2,3);                  % J * (cov_pars * J')
        end

        function periods = get.periods(obj)
            periods = repmat(obj.const_to_periods(), [obj.ncomponents, 1]);
        end

        function omega = get_omega(obj)
            omega = 2*pi./obj.periods;
        end

        function val = get_npars_tid(obj)
            val = ones(1, obj.get_ncomponents).*(1 + 2*obj.nconstituents);
        end

    end
    methods(Access=protected)
        function mustBeValidComponent(obj, components)
            all_comps = obj.get_component_names;
            assert(all(ismember(components, all_comps)),...
                ['Components should be one or more of ''',...
                strjoin(all_comps),...
                ''''])
        end
        function mustBeValidConstituent(obj, constituents)
            all_const = [{'M0'}, obj.supported_constituents];
            assert(all(ismember(constituents, all_const)),...
                'Not all given constituents are supported')
        end
        function mustBeSinOrCos(~, parameter)
            all_pars = {'cos', 'sin'};
            assert(all(ismember(parameter, all_pars)),...
                'Value must be ''sin'', or ''cos''')
        end

        function names = get_names(obj)

            % Forms a cell array of dimensions 1xobj.ncomponents
            % Elements of the cell array are cell arrays containing the
            % tidal expansion names per component

            const_names = obj.constituents; % cell array of length n_comp containing constituents of each component
            names = cell([obj.ncomponents, 1]);
            for comp = 1:obj.ncomponents
                names{comp}{1} = [obj.component_names{comp}, ': M0']; % Subtidal part
                idx = 2;
                for const = 1:length(const_names)
                    names{comp}{idx} = [obj.component_names{comp}, ': ', const_names{const}, 'a']; % Cosine part
                    names{comp}{idx+1} = [obj.component_names{comp}, ': ', const_names{const}, 'b']; % Cosine part
                    idx = idx + 2;
                end
            end
        end

        function val = get_npars(obj)
            val = get_npars_tid(obj);
        end

        function val = get_nconstituents(obj)
            val = numel(obj.constituents);
        end

        function periods = const_to_periods(obj)
            % Returns vector of tidal constituent periods in seconds.
            periods = zeros([1, length(obj.constituents)]);
            for const = 1:numel(obj.constituents)
                switch obj.constituents{const}
                    case 'M2'
                        T = 12.4206012;
                    case 'S2'
                        T = 12;
                    case 'N2'
                        T = 12.65834751;
                    case 'K1'
                        T = 23.93447213;
                    case 'M4'
                        T = 6.210300601;
                    case 'O1'
                        T = 25.81933871;
                    case 'M6'
                        T = 4.140200401;
                    case 'MK3'
                        T = 8.177140247;
                    case 'S4'
                        T = 6;
                    case 'MN4'
                        T = 6.269173724;
                    otherwise
                        error('Unknown tidal constituent')
                end
                periods(1,const) = T*3600; % in seconds
            end
        end
    end

end
