classdef HorizontalPositionFromBottomTracking < ADCPHorizontalPosition
    properties
        description = 'Bottom tracking';
        initial_position (2,1) double {mustBeReal} = [0;0];
        interpolate_missing (1,1) logical = true;
    end
    methods (Access=protected)
        function pos = get_horizontal_position(obj, adcp)
            arguments(Input)
                obj  (1,1) HorizontalPositionFromBottomTracking
                adcp (:,:) VMADCP
            end
            arguments(Output)
                pos (2,:) double
            end
            svel_provider = ShipVelocityFromBT;
            shipvel = svel_provider.ship_velocity(adcp,...
                CoordinateSystem.Earth);
            shipvel = permute(shipvel(:, :, 1 : 2), [3, 2, 1]);
            fnan = any(isnan(shipvel),1);
            if obj.interpolate_missing
                shipvel(:, fnan) =...
                    interp1( ...
                    adcp.time(~fnan)',...
                    shipvel(:,~fnan)',...
                    adcp.time(fnan)')';
            else
                shipvel(:, fnan) = repmat([0; 0],[1, sum(fnan,2)]);
            end
            sv_avg = .5 * (shipvel(:, 1 : end - 1) + shipvel(:, 2 : end) );
            dt = seconds(diff(adcp.time));
            sv_dpos = sv_avg .* dt;
            pos = cumsum([obj.initial_position, sv_dpos],2);
        end
        function tf=get_has_data(~, adcp)
            tf = ismethod(adcp,'btvel');
        end
    end
end