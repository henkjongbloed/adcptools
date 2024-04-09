classdef DataModel <...
        helpers.ArraySupport
    % Model for Cartesian velocity components
    %
    %   This base class implements a dummy model, i.e. each velocity is
    %   represented by the mean. Subclasses can implement linear models for the
    %   velocity. They should implement a linear mapping between model
    %   parameters and cartesian velocity components. Subclasses can
    %   reimplement the get_model methods and the get_npars methods to
    %   implement a new velocity model.
    %
    %   DataModel properties:
    %   ncomponents - number of components in the model
    %   npars - number of parameters for each component
    %   names - names of parameters per component
    %   all_names - all model parameter names
    %
    %   DataModel methods:
    %   get_data - reconstructs data from model parameters
    %
    %  see also: VelocityModel, TaylorModel, TidalModel, DataSolver

    properties(Dependent, SetAccess=private)
        % DataModel/npars (read only) number of parameters in the model
        %
        %   1xNCOMP row vector holding the number of model parameters for each
        %   data components. 
        %
        %   see also: DataModel, get_model
        npars

        % DataModel/ncomponents (read only) components in model
        %
        %   scalar value holding the number of components considered in the
        %   model
        %
        %   see also: DataModel, get_model
        ncomponents

        % DataModel/component_names
        %
        %   ncomponent x 1 cell holding the name of components in the model.
        %
        %   see also: DataModel, names
        component_names
        
        % DataModel/names (read only) of parameters per component
        %
        %   ncomponents x 1 cell holding the names of the parameters in the
        %   model for each component.
        %
        % see also: DataModel, all_names
        names

        % DataModel/all_names (read only) of parameters
        %
        %   1 x sum(npars) cell array with the names of all the parameters
        %   in the model
        %
        % see also: DataModel, names
        all_names
    end
    methods
        function [dat, cov_dat] = get_data(obj, parameters, options)
            %pars, cov_pars,...
             %   n_bvels, d_time, d_s, d_n, d_z, d_sigma)
            % Compute data values from model parameters.
            %
            %   dat = get_data(obj, pars) computes data based on model
            %   parameters.
            %
            %   [dat, cov_dat] = get_data(obj, pars, cov_pars) also compute
            %   covariance matrix of data based on covariance of model
            %   parameters.
            %
            %   [dat, cov_dat] = get_velocity(obj, pars, cov_pars, d_time, d_s,
            %       d_n, d_z, d_sigma) optionally specify the coordinate
            %       offsets at which data will be computed. These default
            %       to zeros if not given. They are all column vectors with
            %       same number of rows as the rows in pars.
            %
            %   see also: DataModel, Solver
            arguments
                obj
                parameters (:,:,:) double
                options.covariance (:,:,:,:) double =...
                    zeros(size(parameters,[1, 2, 2, 3]));
                options.delta_t (:,1) datetime = ...
                    datetime(zeros(size(parameters,1),1), 'ConvertFrom', 'datenum');
                options.delta_s (:,1) double = zeros(size(parameters,1),1);
                options.delta_n (:,1) double = zeros(size(parameters,1),1);
                options.delta_z (:,1) double = zeros(size(parameters,1),1);
                options.delta_sig (:,1) double =...
                    zeros(size(parameters,1),1);
                
            end
            
            M = obj.get_model(...
                options.delta_t,...
                options.delta_s,...
                options.delta_n,...
                options.delta_z,...
                options.delta_sig...
                );

            % reconstruct model matrix with kron product
            np = obj.npars;
            nin = size(parameters,1);
            nc = obj.ncomponents;
            Mnew = zeros(nin,nc,sum(np));
            cum_pars = cumsum([0 np]);
            for c_comp = 1:obj.ncomponents
                Mnew(:,c_comp,cum_pars(c_comp)+1:cum_pars(c_comp+1))=...
                    M(:,1:np(c_comp),c_comp);
            end
            M=Mnew;

            % apply model to obtain vel
            pars = permute(parameters,[1 3 2 4]);
            M = permute(M,[1 4 2 3]);
            dat = helpers.matmult(M,pars);
            dat = permute(dat,[1 3 2]);

            % apply model to obtain covariance matrix
            if nargout > 1
                cov_pars = permute(options.covariance, [1 4 2 3]);
                cov_dat = helpers.matmult(cov_pars, permute(M,[1,2,4,3]));
                cov_dat = helpers.matmult(M, cov_dat);
                cov_dat = permute(cov_dat,[1 3 4 2]);
            end
        end


        function [dat, cov_dat] = get_data2(obj, parameters, options)
            %   dat = get_data(obj, pars) computes data based on model
            %   parameters.
            %
            %   [dat, cov_dat] = get_data(obj, pars, cov_pars) also compute
            %   covariance matrix of data based on covariance of model
            %   parameters.
            %
            %   [dat, cov_dat] = get_velocity(obj, pars, cov_pars, d_time, d_s,
            %       d_n, d_z, d_sigma) optionally specify the coordinate
            %       offsets at which data will be computed. These default
            %       to zeros if not given. They are all column vectors with
            %       same number of rows as the rows in pars.
            %
            %   see also: DataModel, Solver
            arguments
                obj
                parameters (:,:,:) double
                options.covariance (:,:,:,:) double =...
                    zeros(size(parameters,[1, 2, 2, 3]));
                options.delta_t (:,1) datetime = ...
                    datetime(zeros(size(parameters,1),1), 'ConvertFrom', 'datenum');
                options.delta_s (:,1) double = zeros(size(parameters,1),1);
                options.delta_n (:,1) double = zeros(size(parameters,1),1);
                options.delta_z (:,1) double = zeros(size(parameters,1),1);
                options.delta_sig (:,1) double =...
                    zeros(size(parameters,1),1);
                
            end
            
            M = obj.get_model(...
                options.delta_t,...
                options.delta_s,...
                options.delta_n,...
                options.delta_z,...
                options.delta_sig...
                );

            % reconstruct model matrix with kron product
            np = obj.npars;
            nin = size(parameters,1);
            nc = obj.ncomponents;
            Mnew = zeros(nin,nc,sum(np));
            cum_pars = cumsum([0 np]);
            for c_comp = 1:obj.ncomponents
                Mnew(:,c_comp,cum_pars(c_comp)+1:cum_pars(c_comp+1))=...
                    M(:,1:np(c_comp),c_comp);
            end
            M=Mnew;

            % apply model to obtain vel
            pars = permute(parameters,[1 3 2 4]);
            M = permute(M,[1 4 2 3]);
            dat = helpers.matmult(M,pars);
            dat = permute(dat,[1 3 2]);

            % apply model to obtain covariance matrix
            if nargout > 1
                cov_pars = permute(options.covariance, [1 4 2 3]);
                cov_dat = helpers.matmult(cov_pars, permute(M,[1,2,4,3]));
                cov_dat = helpers.matmult(M, cov_dat);
                cov_dat = permute(cov_dat,[1 3 4 2]);
            end
        end


        function val=get.npars(obj)
            val=obj.get_npars();
        end
        function val=get.ncomponents(obj)
            val=obj.get_ncomponents();
        end
        function names = get.names(obj)
            names = obj.get_names;
            assert(iscell(names),...
                'Parameter names should be given as a cell')
            assert(numel(names) == obj.ncomponents, ...
                ['Number of cell elements should match number of',...
                ' components'])
            assert(all(cellfun(@(x) iscellstr(x) || isstring(x), names))...
                , ['Names for each component should be given as cell',...
                ' of character vectors or as strings'])
            assert(isequal(...
                reshape(cellfun(@numel, names),1,[]),...
                obj.npars),...
                ['Number of parameter names for each component should',...
                ' match the number of parameters for each component.'])
        end
        function val = get.component_names(obj)
            val = obj.get_component_names();
            assert(iscellstr(val), ...
                'components must be given as a cell of chars') %#ok<ISCLSTR>
            assert(numel(val)==obj.ncomponents,...
                ['number of component names should match number of ',...
                'components'])
        end
        function val = get.all_names(obj)
            val = [obj.names{:}];
        end
    end
    methods(Abstract)
            % build velocity model matrices
            %
            %   M is NEnsembles x NPars x Ncomp
            %
            %   [Mu, Mv, Mw] = get_model(obj, d_time, d_s, d_n, d_z,
            %   d_sigma) given the offsets of the velocity data to the mesh
            %   cell coordinates return the model matrices which describe the
            %   relation between model parameters and cartesian velocity.
            %   For each cartesian velocity component a model is given in Mu,
            %   Mv and Mw for the x, y, and z component respectively. Mu, Mv,
            %   and Mw are matrices with each row corresponding to a measured
            %   velocity component and each column being the known coefficient
            %   for a model parameter.
            %   For DataModel, the matrices hold a column of ones, meaning
            %   that each measured velocity component is an estimate of the
            %   final mean velocity.
            %   Reimplementing this in a subclass allows to define other, more
            %   complex, but linear velocity models.
            %
            %   d_time, d_s, d_n, d_z, d_sigma are the distances of the
            %   measured velocity data from the center of the cell, or the time
            %   difference with the mesh time.
            %
            %   see also: Solver, TaylorModel, TidalModel.
            get_model(obj, d_time, ~, ~, ~, ~)
    end
    methods(Abstract, Access=protected)
        get_names
        get_ncomponents
        get_npars
        get_component_names
    end
end