classdef ADCP < handle &...
        helpers.ArraySupport
% Abstract base class to represent ADCP datasets
%
%   A = ADCP(...) based on the class of the passed arguments
%   configuration properties are set:
%   Filter - filter
%   HorizontalPositionProvider - horizontal_position_provider
%   VerticalPositionProvider - vertical_position_provider
%
%   ADCP properties (configuration):
%   filters - filters to exclude data from profiles
%   timezone - specifies the timezone used
%   horizontal_position_provider - provides horizontal positioning
%   vertical_position_provider - provides vertical positioning
%   tilts_provider - provides pitch and roll angles
%   heading_provider - provides heading angle
%   instrument_matrix_provider - provides instrument to beam matrix
%   transducer - transducer acoustic properties
%   water - water acoustic properties
%
%   ADCP data properties (read only):
%   nbeams - number of acoustic beams
%   nensembles - number of ensembles
%   ncells - number of depth cells
%   coordinate_system - coordinate system of velocity data
%   beam_angle - angle of acoustic beams with vertical (degrees)
%   pitch - pitch angle (degrees)
%   roll - roll angle (degrees)
%   heading - heading angle (degrees)
%   cellsize - depth cell size (m)
%   time - time and date of ensembles
%   distmidfirstcell - distance to center of first depth cell (m)
%   depth_cell_slant_range - distance along ac. beam to depth cells (m)
%   temperature - instrument temperature (Celsius)
%   salinity - salinity of water (psu)
%   pressure - pressure of water (Pa)
%   echo - received echo intensity (dB)
%   backscatter - volume backscatter strength (dB)
%   horizontal_position - horizontal position of the instrument (m)
%   vertical_position - vertical position of the instrument (m)
%   beam_2_instrument_matrix - beam to instrument transformation matrix
%   instrument_2_beam_matrix - instrument to beam transformation matrix
%
%   ADCP methods:
%   bad - return a mask with ones for data marked as being bad
%   plot_orientations - plots pitch, roll and heading of the ADCP
%   plot_velocity - plot the velocity in Earth coordinates
%   plot_backscatter - plot the volume backscatter strength
%   plot_all - plot all available plots
%   depth_cell_offset - get east,north,up offset of depth cells from ADCP
%   depth_cell_position - get east,north,up position of depth cells
%   
%   ADCP methods (abstract):
%   xform - get transformation matrices
%   velocity - get velocity

    properties
        % ADCP/description
        %
        %   brief description of the ADCP data. Value must be a string
        %   scalar
        %
        %   see also: ADCP
        description(1,1) string = "";

        % ADCP/filters
        %
        % Filters to exclude data from profile data. The filters are given
        % as a vector of object that are subclasses of Filter
        %
        % see also: ADCP, Filter
        filters(:,1) Filter {mustBeNonempty}=Filter;
        
        % ADCP/timezone
        %
        %   Specify the timezone of the data as char. Default is '', i.e.
        %   unzoned. Examples are 'UTC', 'Europe/Amsterdam'
        %
        % see also: ADCP, datetime/TimeZone, timezones
        timezone(1,:) char = ''

        % ADCP/horizontal_position_provider
        %
        %   ADCPHorizontalPosition object specifying the position of the ADCP.
        %
        % see also: ADCP
        horizontal_position_provider(1,:) ADCPHorizontalPosition = ADCPFixedHorizontalPosition;
        
        % ADCP/vertical_position_provider
        %
        %   ADCPVerticalPosition object specifying the position of the ADCP.
        %
        % see also: ADCP
        vertical_position_provider(1,1) ADCPVerticalPosition = ADCPFixedVerticalPosition;

        % ADCP/tilts_provider
        %
        % TiltsProvider object which returns the heading and tilts of the 
        % ADCP.
        %
        % see also: ADCP
        tilts_provider(:,1) TiltsProvider = rdi.TiltsInternal;

        % ADCP/heading_provider
        %
        % HeadingProvider object which returns the heading and tilts of the 
        % ADCP.
        %
        % see also: ADCP
        heading_provider(:,1) HeadingProvider = rdi.HeadingInternal


        % ADCP/instrument_matrix_provider
        %
        %   Specifies the sources for the transformation matrix as a
        %   InstrumentMatrixProvider. This is a vector which allows to
        %   define different sources. The first object capable of providing
        %   the matrix will be used.
        %
        % see also: ADCP, InstrumentMatrixProvider
        instrument_matrix_provider (:,1) InstrumentMatrixProvider =...
            rdi.InstrumentMatrixFromBAngle;

        % ADCP/water_level_object
        %
        %   Specifies the source for the water level. Should be a scalar
        %   WaterLevel object. It defaults to ConstantWaterLevel with a
        %   level set to 0; This property is used for the computation of
        %   the water_level property.
        %
        % see also: ADCP, ADCP.water_level
        water_level_object(1,1) WaterLevel = ConstantWaterLevel(0)
    end
    properties(Access=protected)
        override_transducer
        override_water
    end
    properties(Dependent)
        % ADCP/tranducer
        %
        % returns an acoustics.PistonTransducer object. This object can be
        % modified, and will be used in computations. To reset object to
        % its initial values use the reset_transducer method. Setting the
        % raw property, or the water property will reset the tranducer
        % property
        %
        % see also: ADCP, acoustics.PistonTransducer
        transducer(:,1) acoustics.PistonTransducer 
        
        % ADCP/water
        %
        %   acoustics.Water object specifying the Water characteristics.
        %   use reset_water method to reset the water property to the
        %   values read in the raw adcp data. Changing the raw property
        %   will also reset the water property.
        %
        % see also: ADCP, acoustics.Water
        water(1,1) acoustics.Water;    
    end
    properties(Dependent, SetAccess = private)
        % Number of beams in use by the instrument.
        %
        % see also: ADCP
        nbeams
        
        % Number of ensembles.
        %
        % see also: ADCP
        nensembles
        
        % Number of depth cells.
        %
        % see also: ADCP
        ncells
        
        % Coordinate system used. Returns CoordinateSystems objects.
        %
        % see also: ADCP, CoordinateSystem
        coordinate_system
        
        % Angle of acoustic beam makes with vertical in degrees.
        %
        % see also: ADCP
        beam_angle
        
        % Pitch angle in degrees.
        %
        % see also: ADCP
        pitch
        
        % Roll angle in degrees.
        %
        % see also: ADCP
        roll
        
        % Heading angle in degrees.
        %
        % see also: ADCP
        heading

        % Vertical size of the depth cells (m).
        %
        % see also: ADCP
        cellsize

        % Time the ensembles were measured. Returns datetime objects.
        %
        % see also: ADCP, datetime
        time
        
        % Vertical distance to the middle of the first cell (m).
        %
        % see also: ADCP
        distmidfirstcell
        
        % Slant range, i.e. distance along acoustic beam to depth cells
        %
        % see also: ADCP
        depth_cell_slant_range
        
        % Temperature of the instrument (Celsius)
        %
        % see also: ADCP
        temperature
        
        %  salinity (psu)
        %
        % see also: ADCP
        salinity
        
        % pressure (Pa)
        %
        % see also: ADCP
        pressure

        % Received raw echo intensity (dB).
        %
        % see also: ADCP
        echo
        
        % Received Volume Backscatter strength (dB) computed according to
        % Deines (1999) with the corrections of Gostiaux and van Haren and 
        % the correction in the FSA-031. The backscatter strength is not 
        % corrected for attenuation due to sediment.
        %
        % see also: ADCP, acoustics, Sv2SSC
        backscatter
        
        %   2xN matrix holding the x and y coordinates of the
        %   ADCP in m.
        %
        %   see also: adcp, vertical_position, horizontal_position_provider
        horizontal_position
        
        %   1xN vector holding the z coordinates of the ADCP in m.
        %
        %   see also: adcp, horizontal_position, vertical_position_provider
        vertical_position

        % matrix transforming beam coordinates to instrument coordinates
        %
        %   see also: ADCP, instrument_matrix_provider
        beam_2_instrument_matrix

        % matrix transforming instrument coordinates to beam coordinates
        %
        %   see also: ADCP, instrument_matrix_provider
        instrument_2_beam_matrix

        % matrix holding unit vectors pointing in direction of beams
        %
        %   size: 1 x nensembles x nbeams x 3 components
        %   for each beam it holds the 3 components of the unit vector
        %   pointing along the beam and away from the ADCP. Note that this
        %   matrix can be defined for an uplooking or downlooking position
        %   depending on the conventions of the manufacturer, i.e. the
        %   position for which all tilts are zero.
        %
        %   see also: ADCP, instrument_matrix_provider
        beam_orientation_matrix

        % vertical distance to cells for an untilted ADCP
        %
        %   size: ncells x nensembles x nbeams
        %   holds the vertical position of the depth cells for an
        %   untilted ADCP. If an untilted ADCP is downlooking according to
        %   the manufacturer's conventions this will hold negative values.
        %
        %   see also: ADCP, beam_orientation_matrix
        vertical_range_to_cell

        % water_level
        %
        %   size: 1 x nensembles
        %   returns the water level. The way the water level is computed
        %   can be modified by changing the water_level_object property.
        %
        %   see also: ADCP, water_level_object
        water_level

        % depth_cell_position
        %
        %   size: ncells x nensembles x nbeams x 3
        %   position vector in Earth coordinates
        %
        % see also: ADCP, position, depth_cell_offset
        depth_cell_position
    end
    methods
        %%% Constuctor
        function obj = ADCP(varargin)
            obj.filters=Filter;
            for ca=1:nargin            
                if isa(varargin{ca},'Filter')
                    obj.filters=[obj.filters; varargin{ca}];
                elseif isa(varargin{ca},'ADCPVerticalPosition')
                    obj.vertical_position_provider=varargin{ca};
                elseif isa(varargin{ca},'ADCPHorizontalPosition')
                    obj.horizontal_position_provider=varargin{ca};
                end
            end
        end
        
        %%% Set and get methods
        function val = get.water(obj)
            if isempty(obj.override_water)
                val=acoustics.Water;
                val.temperature=obj.temperature;
                val.salinity=obj.salinity*1000;
            else
                val=obj.override_water;
            end
        end
        function wl=get.water_level(obj)
            wl=obj.water_level_object.get_water_level(obj.time);
        end
        function set.water(obj, val)
            if isempty(val)
                obj.water = acoustics.Water.empty;
                return
            end
            assert(isa(val, 'acoustics.Water'), 'Must be an acoustics.Water object')
            assert(isscalar(val), 'Value must be scalar')
            obj.override_water = val;
        end
        function val = get.transducer(obj)
            if isempty(obj.override_transducer)
                val = obj.get_transducer;
            else
                val = obj.override_transducer;
            end
        end
        function set.transducer(obj, val)
            if isempty(val)
                obj.transducer = acoustics.PistonTransducer.empty;
                return
            end
            assert(isa(val, 'acoustics.Transducer'), 'Must be an acoustics.Transducer object')
            assert(isscalar(val), 'Value must be scalar')
            obj.override_transducer = val;
        end
        function val = get.nbeams(obj)
            val = obj.get_nbeams;
        end
        function val = get.nensembles(obj)
            val = obj.get_nensembles;
        end
        function val = get.ncells(obj)
            val = obj.get_ncells;
        end
        function val = get.coordinate_system(obj)
            val = obj.get_coordinate_system;
        end
        function val = get.beam_angle(obj)
            val = obj.get_beam_angle;
        end
        function val = get.beam_orientation_matrix(obj)
            val = obj.instrument_matrix_provider.beam_orientation_matrix(obj);
        end
        function val = get.pitch(obj)
            val = obj.tilts_provider.pitch(obj);
        end
        function val = get.roll(obj)
            val = obj.tilts_provider.roll(obj);
        end
        function val = get.heading(obj)
            val = obj.heading_provider.heading(obj);
        end
        function val = get.cellsize(obj)
            val = obj.get_cellsize;
        end
        function val = get.time(obj)
            val = obj.get_time;
        end
        function val = get.distmidfirstcell(obj)
            val = obj.get_distmidfirstcell;
        end
        function rng=get.depth_cell_slant_range(obj)
            beam_m = obj.instrument_matrix_provider.beam_orientation_matrix(obj);
            beam_m(:,:,:,1:2)=[];
            rng=obj.vertical_range_to_cell./beam_m;
        end
        function val = get.temperature(obj)
            val = obj.get_temperature;
        end
        function val = get.salinity(obj)
            val = obj.get_salinity;
        end
        function val = get.pressure(obj)
            val = obj.get_pressure;
        end
        function val = get.echo(obj)
            val = obj.get_echo;
        end
        function val = get.backscatter(obj)
            val = obj.get_backscatter;
        end
        function val=get.horizontal_position(obj)
            val=obj.horizontal_position_provider.horizontal_position(obj);
        end
        function val=get.vertical_position(obj)
            val=obj.vertical_position_provider.get_vertical_position(obj);
        end
        function val = get.beam_2_instrument_matrix(obj)
            val = obj.instrument_matrix_provider.b2i_matrix(obj);
        end
        function val = get.instrument_2_beam_matrix(obj)
            val = obj.instrument_matrix_provider.i2b_matrix(obj);
        end
        function val = get.vertical_range_to_cell(obj)
            val = obj.instrument_matrix_provider.vertical_range_to_cell(obj);
        end
        function bad=bad(obj,filter)
            % Mask for profiled data
            %
            %   isbad=bad(obj) returns a mask with ones for data data is marked bad
            %   by the filters.
            %
            %   isbad=bad(obj,filters) allows to specify custom filters, other than
            %   the ones given in the ADCP object.
            %
            %   see also: ADCP, Filter
            if nargin<2
                filter=obj.filters;
            end
            bad=filter(1).bad(obj);
            for co=2:numel(filter)
                bad=bad | filter(co).bad(obj);
            end
        end
        function varargout = plot_orientations(obj)
            % Plots orientations of instrument
            %
            %   obj.plot_orientations plots the pitch, roll and heading.
            %
            %   [ahp, ahr, ahh] = obj.plot_orientation returns the axes
            %   handles for the pitch, roll adn heading axes respectively.
            %
            % see also: ADCP, plot_velocity, plot_filters, plot_all
            no = numel(obj);
            [axpitch, axroll, axhead] = deal(nan(no,1));
            t=tiledlayout(3, no,'TileSpacing','tight','Padding','compact');
            for ca = 1:no
                axpitch(ca)=nexttile;
                plot(obj(ca).time, obj(ca).pitch)
                title(obj(ca).description);
                set(gca,'XTickLabel',[],'XTickMode','manual')
            end
            for ca = 1:no
                axroll(ca)=nexttile;
                plot(obj(ca).time, obj(ca).roll)
                set(gca,'XTickLabel',[],'XTickMode','manual')
            end
            for ca = 1:no
                axhead(ca)=nexttile;
                plot(obj(ca).time, obj(ca).heading)
                linkaxes([axpitch(ca), axroll(ca), axhead(ca)],'x')
            end
            xlabel(t,'Time');
            ylabel(axpitch(1),'Pitch (degrees)')
            ylabel(axroll(1),'Roll (degrees)')
            ylabel(axhead(1),'Heading (degrees)')
            if nargout > 0
                varargout{1} = axpitch;
            end
            if nargout > 1
                varargout{2} = axroll;
            end
            if nargout > 2
                varargout{3} = axhead;
            end
        end
        function hfout=plot_velocity(obj,vel)
            % plot the velocity profiles
            %
            %   plot_velocity(obj) plot the velocities in Earth coordinate
            %   system
            %
            %   see also: ADCP
            no = numel(obj);
            t=tiledlayout(3, no,'TileSpacing','tight','Padding','compact');
            vel_pos={obj.depth_cell_position};
            vel_pos=cellfun(@(x) mean(x(:,:,:,3),3,'omitnan'), vel_pos,...
                'UniformOutput', false);
            if nargin < 2
                vel={obj.velocity(CoordinateSystem.Earth)};          
            end
            [axx, axy, axz] = deal(nan(1,no));
            tim = cell(1,no);
            hf=gcf;
            for co = 1:no
                axx(co)=nexttile;
                tim{co}=obj(co).time;
                pcolor(tim{co},vel_pos{co},vel{co}(:,:,1));
                shading flat
                hc=colorbar;
                ylabel(hc,'V_E (m/s)')
                title(axx(co),obj(co).description)
                set(gca,'XTickLabel',[],'XTickMode','manual')
            end
            for co = 1:no
                axy(co)=nexttile;
                pcolor(tim{co},vel_pos{co},vel{co}(:,:,2));
                shading flat
                hc=colorbar;
                ylabel(hc,'V_N (m/s)')
                set(gca,'XTickLabel',[],'XTickMode','manual')
            end
            for co = 1:no
                axz(co)=nexttile;
                pcolor(tim{co},vel_pos{co},vel{co}(:,:,3));
                shading flat
                hc=colorbar;
                ylabel(hc,'V_U (m/s)')
                linkaxes([axx(co), axy(co), axz(co)],'xy')
            end
            ylabel(t,'Vertical position (m)')
            xlabel(t,'Time')
            if nargout>0
                hfout=hf;
            end
        end
        function hfout=plot_backscatter(obj)
            % plot the backscatter profiles
            %
            %   plot_velocity(obj) plot the backscatter profiles
            %   system
            %
            %   see also: ADCP

            no = numel(obj);
            nb = [obj.nbeams];
            nbmax = max(nb);
            t=tiledlayout(nbmax, no,'TileSpacing','tight','Padding','compact');
            sv_pos={obj.depth_cell_position};
            sv_pos=cellfun(@(x) mean(x(:,:,:,3),3,'omitnan'), sv_pos,...
                'UniformOutput', false);
            axh = deal(nan(nbmax,no));
            tim={obj.time};
            sv={obj.backscatter};
            for cb = 1:nbmax
                for co = 1:no
                    axh(cb,co)=nexttile;
                    if nb(co)<cb
                        continue
                    end
                    pcolor(tim{co},sv_pos{co},sv{co}(:,:,cb));
                    shading flat
                    hc=colorbar;
                    ylabel(hc,['Sv Beam ',num2str(cb),' (dB)'])
                    shading flat                      
                    if cb == nbmax
                        title(axh(1,co),obj(co).description)
                        set(axh(end,co),'XTickLabel',[],'XTickMode','manual')
                        linkaxes(axh(:,co),'xy')
                    end
                end
            end
            xlabel(t, 'Time (s)')
            ylabel(t, 'Vertical position(m)')
            if nargout>0
                hfout=gcf;
            end
        end
        function plot_filters(obj, varargin)
            if ~isscalar(obj)
                tpar = tiledlayout("flow");
                for co = 1:numel(obj)
                    t = tiledlayout(tpar,1,1);
                    t.Layout.Tile = co;
                    obj(co).plot_filters(t, varargin{:})
                    title(t, obj(co).description)
                end
                return
            end
            obj.filters.plot(obj, varargin{:});
        end
        function plot_all(obj)
            figure
            obj.plot_orientations;
            figure
            obj.plot_filters;
            figure
            obj.plot_velocity;
            figure
            obj.plot_backscatter;
        end
        function varargout=depth_cell_offset(obj,varargin)
            % Computes the xyz offset to the profiled depth cells
            %
            %   pos=depth_cell_offset(obj) returns the offset vector (ncells x
            %   nensembles x nbeams x 3) in Earth coordinates
            %
            %   pos = depth_cell_offset(obj,dst) specify the coordinate system
            %   in which the offsets should be returned. dst is a
            %   CoordinateSystem object.
            %
            %   see also: ADCP, depth_cell_position
            if ~isscalar(obj)
                varargout = obj.run_method_single(nargout, 'depth_cell_offset', varargin{:});
                return
            end
            if nargin < 2
                varargin{1}=CoordinateSystem.Earth;
            end
            tm=obj.xform(CoordinateSystem.Instrument,...
                varargin{1}, 'Geometry', true);
            tm(:,:,:,4)=[];
            tm(:,:,4,:)=[];
            beam_mat = obj.instrument_matrix_provider.beam_orientation_matrix(obj);
            tm = helpers.matmult(beam_mat,tm, 3, 4);
            varargout{1}=tm.*obj.depth_cell_slant_range;
        end
        function pos=get.depth_cell_position(obj)
            pos=obj.depth_cell_offset + permute([obj.horizontal_position; obj.vertical_position],[3,2,4,1]);
        end
    end
    methods(Abstract, Access = protected)
        val = get_nbeams(obj)
        val = get_nensembles(obj)
        val = get_ncells(obj)
        val = get_coordinate_system(obj)
        val = get_beam_angle(obj)
        val = get_cellsize(obj)
        val = get_time(obj)
        val = get_distmidfirstcell(obj)
        val = get_temperature(obj)
        val = get_salinity(obj)
        val = get_pressure(obj)
        val = get_echo(obj)
        val = get_backscatter(obj)
        val = get_transducer(obj)
    end
    methods
        % Get transformation matrices for coordinate transformations
        %
        %   tm=xform(obj,dst) get the transformation matrices for the
        %   transformation from obj.coordinate_system to the given
        %   destination coordinate system. matrix will be 1 x nensembles x
        %   4 x 4
        %
        %   tm=xform(obj,dst,src) to specify a custom source coordinate
        %   system
        %
        %   tm=xform(...,'ParameterName', parameterValue) allows to
        %   specify the following options:
        %   'UseTilts' - if unset checks the ADCP data in inverse
        %   transformations to know whether to use the tilts. If set to
        %   true will always use the tilts in invers transformations,
        %   if set to false will never use tilts in inverse
        %   transformations.
        %
        %   see also: ADCP
        function varargout = xform(obj,varargin)
            if ~isscalar(obj)
                varargout = obj.run_method_single(nargout, 'xform', ...
                    varargin{:});
                return
            end
            varargout{1} = obj.get_xform(varargin{:});
         end
        % velocity profile data
        %
        %   vel=velocity(obj) returns the profiled velocity in m/s.
        %
        %   vel=velocity(obj,dst) specify destination coordinate system as
        %   CoordinateSystem object.
        %
        %   vel=velocity(obj,dst,filter) specify filter to be used insted
        %   of the ones specified in the current object.
        %
        %   see also: ADCP, CoordinateSystem
        function varargout = velocity(obj,varargin)
            if ~isscalar(obj)
                varargout = obj.run_method_single(nargout, 'velocity', ...
                    varargin{:});
                return
            end
            varargout{:} = obj.get_velocity(varargin{:});
        end
    end
    methods(Static)
        function inv=invert_xform(tm)
            % Inverts transformation matrices
            %
            %   Inverts transformation matrices given in tm. Matrices are
            %   defined in the third and fourth dimension.
            %
            % see also: ADCP, xform
            inv=nan(size(tm));
            for ce=1:size(tm,2)
                inv(1,ce,:,:)=shiftdim(inv(squeeze(tm(1,ce,:,:))),-2);
            end
        end
    end
    methods(Access=protected, Abstract)
        get_velocity(obj,dst,filter)
        get_xform(obj,dst,src,varargin)
    end
end
