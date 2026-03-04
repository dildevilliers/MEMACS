classdef Vect3D
   %VECT3D   Class of vectors in 3-D space
   % Objects of this class are vectors in 3-D space, described in cartesian,
   % cylindrical and spherical coordinates. Vectors can be translated,
   % rotated, scaled, as well as added and subtracted from each other. They
   % can also be defined and redefined in terms of different
   % CoordinateSystem bases. Several plotting functions are provided.
   
    properties (SetAccess = private)
        v(3,:) double = [0; 0; 1]   % Cartesian vector values: [x;y;z] components
        pStart(1,:) Pnt3D = Pnt3D   % Starting points of all elements
       
    end
    
    properties (Dependent = true)
        N(1,1) double  % Number of vector elements in the object

        abs(1,:) double %  Magnitude of the vectors

        pEnd(1,:) Pnt3D  % End points of the vectors
        pMid(1,:) Pnt3D  % Mid points of the vectors

        % Unit vectors at the starting points in standard coordinate systems
        u_x(3,:) double
        u_y(3,:) double
        u_z(3,:) double
        u_rho(3,:) double
        u_ph(3,:) double
        u_r(3,:) double
        u_th(3,:) double

        % Vectors projected onto the different unit vectors in all the coordinate systems
        v_x(1,:) double
        v_y(1,:) double
        v_z(1,:) double
        v_rho(1,:) double
        v_ph(1,:) double
        v_r(1,:) double
        v_th(1,:) double
        
        

        % th % polar angle in radians
        % ph % azimuth angle in radians
        % el % elevation angle in radians
        % r  % distance from origin
        % rho % distance from z-axis
        
    end
    
    methods
        function obj = Vect3D(v,pStart)
            % VECT3D class constructor
            % obj = Vect3D(v,pStart) Can be empty, which return a z-directed unit-vector at the
            % origin.  
            %
            % Inputs
            % - v:  Matrix of vector values in m ([0;0;1]) ([x;y;z] in m)
            % - pStart:  Pnt3D vector of start positions of the vectors (all at [0;0;0])
            %             Can also be a [3 x Nv] vector of points in m
            %
            % Outputs
            % - obj:  Vect3D object
            %
            % Dependencies
            % -
            %
            % Created: 2026-03-04, Dirk de Villiers
            % Updated: 2026-03-04, Dirk de Villiers
            %
            % Tested : Matlab R2025b, Dirk de Villiers
            %  Level : 0
            %   File : 
            %
            % Example
            %   v = Vect3D()
            %   v.plot
            
            if nargin > 0 && ~isempty(v), obj.v = v; end

            if nargin > 1 && ~isempty(pStart)
                if isa(pStart,'Pnt3D')
                    if pStart.N == 1
                        obj.pStart = Pnt3D(ones(1,obj.N).*pStart.x,ones(1,obj.N).*pStart.y,ones(1,obj.N).*pStart.z);
                    else
                        obj.pStart = pStart;
                    end
                else
                    if size(pStart,2) == 1
                        obj.pStart = Pnt3D(ones(1,obj.N).*pStart(1),ones(1,obj.N).*pStart(2),ones(1,obj.N).*pStart(3));
                    else
                        obj.pStart = Pnt3D(pStart);
                    end
                end
            end

            % Expand singe vector if needed
            if obj.N == 1 && obj.pStart.N > 1
                obj.v = repmat(obj.v, 1, obj.pStart.N);
            end

            assert(obj.N == obj.pStart.N,'pStart and v must have the same number of points (or columns)')
            assert(all(obj.abs > eps),'All vector elements must have a finite magnitude')
        end

        function N = get.N(obj)
            N = size(obj.v,2);
        end

        function abs = get.abs(obj)
            abs = vecnorm(obj.v);
        end

        function pEnd = get.pEnd(obj)
            pEnd = obj.pStart.addVect(obj.v);
        end

        function pMid = get.pMid(obj)
            pMid = Pnt3D((obj.pStart.x+obj.pEnd.x)./2,(obj.pStart.y+obj.pEnd.y)./2,(obj.pStart.z+obj.pEnd.z)./2);
        end

        function u_x = get.u_x(obj)
            u_x = repmat([1;0;0],1,obj.N);
        end

        function u_y = get.u_y(obj)
            u_y = repmat([0;1;0],1,obj.N);
        end

        function u_z = get.u_z(obj)
            u_z = repmat([0;0;1],1,obj.N);
        end

        function u_rho = get.u_rho(obj)
            u_rho = [cos(obj.pStart.ph);...
                     sin(obj.pStart.ph); zeros(1,obj.N)];   
        end

        function u_ph = get.u_ph(obj)
            u_ph = [-sin(obj.pStart.ph);...
                     cos(obj.pStart.ph); zeros(1,obj.N)];   
        end

        function u_r = get.u_r(obj)
            u_r = [sin(obj.pStart.th).*cos(obj.pStart.ph);...
                   sin(obj.pStart.th).*sin(obj.pStart.ph);...
                   cos(obj.pStart.th)];   
        end

        function u_th = get.u_th(obj)
            u_th = [cos(obj.pStart.th).*cos(obj.pStart.ph);...
                    cos(obj.pStart.th).*sin(obj.pStart.ph);...
                    -sin(obj.pStart.th)];   
        end

        function v_x = get.v_x(obj)
            v_x = sum(obj.v.*obj.u_x);
        end

        function v_y = get.v_y(obj)
            v_y = sum(obj.v.*obj.u_y);
        end

        function v_z = get.v_z(obj)
            v_z = sum(obj.v.*obj.u_z);
        end

        function v_rho = get.v_rho(obj)
            v_rho = sum(obj.v.*obj.u_rho);
        end

        function v_ph = get.v_ph(obj)
            v_ph = sum(obj.v.*obj.u_ph);
        end

        function v_r = get.v_r(obj)
            v_r = sum(obj.v.*obj.u_r);
        end

        function v_th = get.v_th(obj)
            v_th = sum(obj.v.*obj.u_th);
        end

       
        
        %% Property setters
        
         %% Object operations
        function obj = getNvects(obj,I)
            %GETvects   Returns a subset of the vectors.
            % obj = getNvects(obj,I) returns a Vect3D object of only the vectors 
            % in the requested indexes
            %
            % Inputs
            % - obj:  Vect3D object
            % - I:    Int array of required indexes 
            %
            % Outputs
            % - obj:  Vect3D object
            %
            % Dependencies
            % - 
            %
            % Created: 2026-03-04, Dirk de Villiers
            % Updated: 2026-03-04, Dirk de Villiers
            %
            % Tested : Matlab R2025b, Dirk de Villiers
            %  Level : 0
            %   File : 
            %
            
            assert(isscalar(obj),'Can only get N points for a single object - not a vector of objects')
            obj = Vect3D(obj.v(:,I),obj.pStart(:,I));
        end
        
        %% Overloaded methods
       

        %% Object element operations    
        function obj = scale(obj,scaleVal)
            %SCALE   Scale the vector lengths
            % obj = scale(obj,scaleVal) scales the vectors of
            % the object by scaleVal, but keeps the starting points unchanged
            %
            % Inputs
            % - obj:  Vect3D object
            % - scaleVal: Scaling factor. Can be scalar or numel(scaleVal) == 3
            %
            % Outputs
            % - obj:  Vect3D object
            %
            % Dependencies
            % - 
            %
            % Created: 2026-03-04, Dirk de Villiers
            % Updated: 2026-03-04, Dirk de Villiers
            %
            % Tested : Matlab R2025b, Dirk de Villiers
            %  Level : 0 
            %   File : 
            %
            
            assert(isscalar(scaleVal),'scaleVal must be a finite scalar');
            
            obj.v = obj.v.*scaleVal;
        end
        
        
        function obj = changeBase(obj,coor_new,coor_base)
            %CHANGEBASE   Changes the base coordinate system.
            % obj = changeBase(obj,coor_new,coor_base) changes the base
            % coordinate system in which the object points and vectors where defined
            % from coor_base to coor_new. If only 2 parameters are
            % provided, coor_base is assumed to be the global
            % coordinate system. This method is often used to get points/vectors in
            % the global coordinate system that were defined in some other
            % local system. 
            %
            % Inputs
            % - obj:        Vect3D object
            % - coor_new:   New CoordinateSystem object 
            % - coor_base:  Original/Base CoordinateSystem object 
            %
            % Outputs
            % - obj:  Vect3D object
            %
            % Dependencies
            % - 
            %
            % Created: 2026-03-04, Dirk de Villiers
            % Updated: 2026-03-04, Dirk de Villiers
            %
            % Tested : Matlab R2025b, Dirk de Villiers
            %  Level : 2 
            %   File : testScript_Vect3D.m
            %
            
            
            if nargin < 3
                coor_base = [];
            end
            
            % Split up like this for speed if one of the coordinate inputs are empty (which happens a lot)
            if isempty(coor_base)
                Vglob = obj.v;
            else
                Vbase = obj.v;
                % Rotate the vectors to the global reference
                Q = dirCosine([],coor_base);
                Vglob = Q\Vbase;
            end
            % Go from global to coorNew
            if isempty(coor_new)
                Vprime = Vglob;
            else
                % Rotate the vectors to the new reference
                Q = dirCosine(coor_new,[]);
                Vprime = Q\Vglob;
            end
            
            % Make the object
            obj.v = Vprime;
            obj.pStart = obj.pStart.changeBase(coor_new,coor_base);
        end

        %% Plotting
        function plot(obj,varargin)
            %PLOT Plots the vectors at their starting points.
            %
            % Inputs
            % - obj: Vect3D object
            % * Arbitrary number of pairs of arguments: ...,keyword,value,... where
            %   keywords and values are from the sets
            %   -- lineStyle:  Standard MATLAB linestyle chars ('-')
            %   -- lineColor: Standard MATLAB color representations ('k')
            %   -- lineWidth:  Positive real number (1)
            %
            % Outputs
            % - []
            %
            % Dependencies
            % - 
            %
            % Created: 2026-03-04, Dirk de Villiers
            % Updated: 2026-03-04, Dirk de Villiers
            %
            % Tested : Matlab R2025b, Dirk de Villiers
            %  Level : 0
            %   File : 
            %
            
            parseobj = inputParser;
            parseobj.FunctionName = 'plot';
            
            typeValidationObj = @(x) validateattributes(x,{'Vect3D'},{'numel',1},'plot','obj',1);
            addRequired(parseobj,'obj',typeValidationObj);

            typeValidationLineStyle = @(x) validateattributes(x,{'char'},{},'plot','lineStyle');
            addParameter(parseobj,'lineStyle','-',typeValidationLineStyle);
            
            typeValidationLineColor = @(x) validateattributes(x,{'char','double'},{},'plot','lineColor');
            addParameter(parseobj,'lineColor','k',typeValidationLineColor);
            
            typeValidationLineWidth = @(x) validateattributes(x,{'double'},{'real','positive'},'plot','lineWidth');
            addParameter(parseobj,'lineWidth',1,typeValidationLineWidth);
            
            parse(parseobj, obj, varargin{:});
            
            lineStyle = parseobj.Results.lineStyle;
            lineColor = parseobj.Results.lineColor;
            lineWidth = parseobj.Results.lineWidth;
            
            % Check the input sizes
            x1 = obj.pStart.x;
            y1 = obj.pStart.y;
            z1 = obj.pStart.z;
            
            % x2 = x1 + obj.v(1,:);
            % y2 = y1 + obj.v(2,:);
            % z2 = z1 + obj.v(3,:);
            % plot3([x1;x2],[y1;y2],[z1;z2],'linestyle',lineStyle,...
            %     'color',lineColor,'lineWidth',lineWidth);
            quiver3(x1,y1,z1,obj.v(1,:),obj.v(2,:),obj.v(3,:),LineWidth=lineWidth,LineStyle=lineStyle,Color=lineColor)
            grid on, hold on
            xlabel('x (m)')
            ylabel('y (m)')
            zlabel('z (m)')

            obj.pStart.plot('marker','o','markerEdgeColor',lineColor)
            obj.pEnd.plot('marker','.','markerEdgeColor',lineColor)
            
        end
        
        % function plot(obj,varargin)
        %     %PLOT   Plots a Pnt3D object.
        %     % plot(obj,varargin) plots all the points in the object on a
        %     % 3-D graph. There are several options for the plotting style,
        %     % all a subset of the standard MATLAB plot3 options
        %     %
        %     % Inputs
        %     % - obj: Pnt3D object
        %     % * Arbitrary number of pairs of arguments: ...,keyword,value,... where
        %     %   keywords and values are from the sets
        %     %   -- marker:      Standard MATLAB marker chars ('.')
        %     %   -- markerEdgeColor: Standard MATLAB color representations ('k')
        %     %   -- markerFaceColor: Standard MATLAB color representations ('none')
        %     %   -- markerSize:  Positive real number (10)
        %     %   -- lineStyle:  Standard MATLAB linestyle chars ('none')
        %     %   -- lineColor: Standard MATLAB color representations ('k')
        %     %   -- lineWidth:  Positive real number (1)
        %     %
        %     % Outputs
        %     % - []
        %     %
        %     % Created: 2019-05-07, Dirk de Villiers
        %     % Updated: 2020-03-23, Dirk de Villiers
        %     %
        %     % Tested : Matlab R2018b, Dirk de Villiers
        %     %  Level : 2
        %     %   File : testScript_Pnt3D
        %     %
        %     % Example
        %     %   p = Pnt3D;
        %     %   p.plot
        %     %
        %     %   [x,y] = deal(1:5);
        %     %   z = 0; 
        %     %   p = Pnt3D(x,y,z);
        %     %   p.plot('marker','o','lineStyle','-','lineWidth',2)
        % 
        %     parseobj = inputParser;
        %     parseobj.FunctionName = 'plot';
        % 
        %     typeValidationObj = @(x) validateattributes(x,{'Pnt3D'},{},'plot','obj',1);
        %     addRequired(parseobj,'obj',typeValidationObj);
        % 
        %     typeValidationMarker = @(x) validateattributes(x,{'char'},{},'plot','marker');
        %     addParameter(parseobj,'marker','.',typeValidationMarker);
        % 
        %     typeValidationMarkerEdgeColor = @(x) validateattributes(x,{'char','double'},{},'plot','markerEdgeColor');
        %     addParameter(parseobj,'markerEdgeColor','k',typeValidationMarkerEdgeColor);
        % 
        %     typeValidationMarkerFaceColor = @(x) validateattributes(x,{'char','double'},{},'plot','markerFaceColor');
        %     addParameter(parseobj,'markerFaceColor','none',typeValidationMarkerFaceColor);
        % 
        %     typeValidationMarkerSize = @(x) validateattributes(x,{'double'},{'real','positive'},'plot','markerSize');
        %     addParameter(parseobj,'markerSize',10,typeValidationMarkerSize);
        % 
        %     typeValidationLineStyle = @(x) validateattributes(x,{'char'},{},'plot','lineStyle');
        %     addParameter(parseobj,'lineStyle','none',typeValidationLineStyle);
        % 
        %     typeValidationLineColor = @(x) validateattributes(x,{'char','double'},{},'plot','lineColor');
        %     addParameter(parseobj,'lineColor','k',typeValidationLineColor);
        % 
        %     typeValidationLineWidth = @(x) validateattributes(x,{'double'},{'real','positive'},'plot','lineWidth');
        %     addParameter(parseobj,'lineWidth',1,typeValidationLineWidth);
        % 
        %     parse(parseobj, obj, varargin{:});
        % 
        %     marker = parseobj.Results.marker;
        %     markerEdgeColor = parseobj.Results.markerEdgeColor;
        %     markerFaceColor = parseobj.Results.markerFaceColor;
        %     markerSize = parseobj.Results.markerSize;
        %     lineStyle = parseobj.Results.lineStyle;
        %     lineColor = parseobj.Results.lineColor;
        %     lineWidth = parseobj.Results.lineWidth;
        % 
        %     if numel(obj) > 1, obj = obj.fuse; end
        % 
        %     plot3(obj.x(:),obj.y(:),obj.z(:),'linestyle',lineStyle,...
        %         'color',lineColor,'lineWidth',lineWidth,...
        %         'marker',marker,'markerEdgeColor',markerEdgeColor,...
        %         'markerFaceColor', markerFaceColor,'markerSize',markerSize), grid on
        %     xlabel('x (m)')
        %     ylabel('y (m)')
        %     zlabel('z (m)')
        % end
        
        %% Splitting and fusing
        function obj = split(obj1)
            %SPLIT splits the single object in a vector of objects
            % obj = split(obj1) takes the input object, which in general
            % has a vector of point values, and splits it into a similar
            % size vector of objects - each with a single element
            %
            % Inputs
            % - obj1:  Pnt3D object
            %
            % Outputs
            % - obj:  Vector of Pnt3D objects
            %
            % Dependencies
            % - 
            %
            % Created: 2020-03-18, Dirk de Villiers
            % Updated: 2020-03-18, Dirk de Villiers
            %
            % Tested : Matlab R2018b, Dirk de Villiers
            %  Level : 1 
            %   File : 
            %
            % Example
            %  x = [1,2,3;4,5,6];
            %  p = Pnt3D(x,0,0)
            %  p1 = p.split
            
            if ~isempty(obj1)
                obj(numel(obj1.x)) = Pnt3D;
                for pp = 1:numel(obj)
                    obj(pp) = Pnt3D(obj1.x(pp),obj1.y(pp),obj1.z(pp));
                end
                obj = reshape(obj,size(obj1.x));
            else
                obj = Pnt3D.empty;
            end
        end
        
        function obj = fuse(obj1)
            %FUSE fuses the matrix of objects into one
            % obj = fuse(obj1) fuses the matrix of single element objects
            % into one object with a matrix of elements
            %
            % Inputs
            % - obj1:  Matrix of Pnt3D objects
            %
            % Outputs
            % - obj:  Pnt3D object with a matrix of elements
            %
            % Dependencies
            % - 
            %
            % Created: 2020-03-18, Dirk de Villiers
            % Updated: 2020-03-18, Dirk de Villiers
            %
            % Tested : Matlab R2018b, Dirk de Villiers
            %  Level : 1 
            %   File : 
            %
            % Example
            %  x = [1,2,3;4,5,6];
            %  for rr = 1:2
            %   for cc = 1:3
            %     p(rr,cc) = Pnt3D(x(rr,cc),0,0);
            %   end
            %  end
            %  p1 = p.fuse
            
            if isempty(obj1)
                obj = obj1;
            else
                xi = reshape([obj1.x],size(obj1));
                yi = reshape([obj1.y],size(obj1));
                zi = reshape([obj1.z],size(obj1));
                assert(numel(xi) == numel(obj1),'Each object in the input array may only have a single element');
                obj = Pnt3D(xi,yi,zi);
            end
        end
        
        %% Deprecated
        function obj = translate(obj,DEL)
            %TRANSLATE   DEPRECATED - USE addVect or plus
            % obj = translate(obj,DEL) translates the [x,y,z] components of
            % the object by DEL, which must be a 3 element
            % vector [xDel,yDel,zDel]. Pretty much the same thing as
            % Pnt3D.plus...
            %
            % Inputs
            % - obj:  Pnt3D object
            % - DEL: Translation vector with numel(DEL)=3 ordered [dx,dy,dz].
            %
            % Outputs
            % - obj:  Pnt3D object
            %
            % Dependencies
            % - 
            %
            % Created: 2019-05-07, Dirk de Villiers
            % Updated: 2019-05-07, Dirk de Villiers
            %
            % Tested : Matlab R2018b, Dirk de Villiers
            %  Level : 1 
            %   File : 
            %
            % Example
            %   [x,y,z] = deal(1:5);
            %   p = Pnt3D(x,y,z);
            %   DEL = [3,3,0.5].';
            %   pd = p.translate(DEL)
            %   p.plot('marker','.'), hold on
            %   pd.plot('marker','o')
            
            warning('Pnt3D.translate is deprecated and not maintained.  Use addVect or plus, depending on the input')
            
            Xp = obj.pointMatrix + DEL(:);
            obj.x = reshape(Xp(1,:),size(obj));
            obj.y = reshape(Xp(2,:),size(obj));
            obj.z = reshape(Xp(3,:),size(obj));
        end
    end
    
    methods (Static = true)
        function obj = sph(PH,TH,R)
            %SPH   Define object in spherical coordinates
            % obj = sph(PH,TH,R) is a constructor for the object in
            % spherical coordinates.  All inputs must be scalar or the same
            % size.
            %
            % Inputs
            % - PH:  Matrix of azimuth angles in radians
            % - TH:  Matrix of polar angles in radians
            % - R:   Matrix of radii in m (default 1)
            %
            % Outputs
            % - obj:  Pnt3D object
            %
            % Dependencies
            % - 
            %
            % Created: 2019-05-07, Dirk de Villiers
            % Updated: 2019-05-07, Dirk de Villiers
            %
            % Tested : Matlab R2018b, Dirk de Villiers
            %  Level : 1 
            %   File : 
            %
            % Example
            %   N = 91;
            %   PH = linspace(0,2*pi,N);
            %   TH = linspace(0,pi/2,N);
            %   R = linspace(1,2,N);
            %   p = Pnt3D.sph(PH,TH,R);
            %   p.plot
            
            assert(all(size(PH)==size(TH)),'PH and TH must be the same size');
            
            if nargin < 3, R = ones(size(PH)); end
            [X,Y,Z] = sph2cart(PH,pi/2 - TH,R);
            obj = Pnt3D(X,Y,Z);
        end
        
        function obj = pol(PH,RHO,Z)
            %POL   Define object in polar coordinates
            % obj = pol(PH,RHO,Z) is a constructor for the object in
            % polar coordinates.  All inputs must be scalar or the same
            % size.  The third input can be omitted for points in the x-y
            % plane
            %
            % Inputs
            % - PH:  Matrix of azimuth angles in radians
            % - RHO: Matrix of distances from the z-axis in m
            % - Z:   Matrix of z-values in m (0)
            %
            % Outputs
            % - obj:  Pnt3D object
            %
            % Dependencies
            % - 
            %
            % Created: 2019-05-07, Dirk de Villiers
            % Updated: 2019-05-07, Dirk de Villiers
            %
            % Tested : Matlab R2018b, Dirk de Villiers
            %  Level : 1 
            %   File : 
            %
            % Example
            %   N = 91;
            %   PH = linspace(0,2*pi,N);
            %   RHO = linspace(0,4,N);
            %   Z = linspace(1,2,N);
            %   p = Pnt3D.sph(PH,RHO,Z);
            %   p.plot
            
            % Define in polar coordinates
            if nargin == 2
                Z = 0;
            end
            [X,Y] = pol2cart(PH,RHO,Z);
            obj = Pnt3D(X,Y,Z);
        end
    
    end
    
    methods (Access = private)
        
    end
    
end
