classdef TaylorTidalScalarModel < TaylorTidalModel & ScalarModel
    methods
        function obj = TaylorTidalScalarModel(varargin)
            obj = obj@TaylorTidalModel(varargin{:});
            obj = obj@ScalarModel(varargin{:});
            order_unass = intersect(obj.unassigned_properties, { ...
                's_order', 'n_order', 'time_order', 'z_order',...
                'sigma_order'});
            for co = 1:numel(order_unass)
                obj.assign_property(order_unass{co}, 0);
            end
        end

        function M = get_model(obj, d_t, d_s, d_n, d_z, d_sigma)
            M = get_model@TaylorTidalModel(obj, d_t, d_s, d_n, d_z, d_sigma);
        end
    end
    methods(Access = protected)
        function val = get_ncomponents(~)
            val = 1;
        end
        function val = get_npars(obj)
            val = get_npars@TaylorTidalModel(obj).*get_npars@ScalarModel(obj);
        end
        function val = get_names(obj)
            val = get_names@TaylorTidalModel(obj);
        end
    end
end