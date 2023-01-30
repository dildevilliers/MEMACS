classdef Pnt3D
   %PNT3D   Class of points in 3-D space
   % Objects of this class are points in 3-D space, described in cartesian,
   % cylindrical and spherical coordinates. Points can be translated,
   % rotated scaled, as well as added and subtracted from each other. They
   % can also be defined and redefined in terms of different
   % CoordinateSystem bases. Several plotting functions are provided.
   
    properties (SetAccess = private)
        x double = 0 % x value in m
        y double = 0 % y value in m
        z double = 0 % z value in m
%         th % polar angle in radians
%         ph % azimuth angle in radians
%         el % elevation angle in radians
%         r  % distance from origin
%         rho % distance from z-axis
    end
    
    properties (Dependent = true)
        th % polar angle in radians
        ph % azimuth angle in radians
        el % elevation angle in radians
        r  % distance from origin
        rho % distance from z-axis
    end
    
    methods
        function obj = Pnt3D(X,Y,Z)
            % PNT3D class constructor
            % obj = Pnt3D(X,Y,Z) Can be empty, which return a point at the
            % origin.  If just one input is provided, it must have 3 rows,
            % which are interpreted as [x;y;z]
            %
            % Inputs
            % - X:  Matrix of x-values in m (0) (or [x;y;z] in m)
            % - Y:  Matrix of y-values in m (0)
            % - Z:  Matrix of z-values in m (0)
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
            %  Level : 2
            %   File : testScript_Pnt3D.m
            %
            % Example
            %   p = Pnt3D([1:3],[2:2:6],4)
            %   p.plot
            
            if nargin == 1
                assert(size(X,1)==3,'One input requires 3 rows')
                shapeVect = size(X);
                if numel(shapeVect) > 2
                    shapeVect = shapeVect(2:end);
                    obj.x = reshape(X(1,:),shapeVect);
                    obj.y = reshape(X(2,:),shapeVect);
                    obj.z = reshape(X(3,:),shapeVect);
                else
                    obj.x = X(1,:);
                    obj.y = X(2,:);
                    obj.z = X(3,:);
                end
            elseif nargin == 3
                % Get all the same size
                obj.x = (Y+eps(realmin))./(Y+eps(realmin)).*(Z+eps(realmin))./(Z+eps(realmin)).*X;
                obj.y = (X+eps(realmin))./(X+eps(realmin)).*(Z+eps(realmin))./(Z+eps(realmin)).*Y;
                obj.z = (X+eps(realmin))./(X+eps(realmin)).*(Y+eps(realmin))./(Y+eps(realmin)).*Z;
            end
        end
        
        function ph = get.ph(obj)
            ph = cart2sph(obj.x,obj.y,obj.z);
        end
        
        function el = get.el(obj)
            [~,el] = cart2sph(obj.x,obj.y,obj.z);
        end
        
        function r = get.r(obj)
            [~,~,r] = cart2sph(obj.x,obj.y,obj.z);
        end
        
        function th = get.th(obj)
            th = pi/2 - obj.el;
        end
        
        function rho = get.rho(obj)
            rho = hypot(obj.x,obj.y);
        end
        
        %% Property setters
        function obj = setX(obj,x)
            %SETX  Set the x-value of the object
            % obj = setX(obj,x)
            assert(numel(obj) == 1,'Can only set a single object - not a vector of objects')
            obj.x = (obj.z+eps(realmin))./(obj.z+eps(realmin)).*(obj.y+eps(realmin))./(obj.y+eps(realmin)).*x;
        end
        
        function obj = setY(obj,y)
            %SETY  Set the y-value of the object
            % obj = setY(obj,y)
            assert(numel(obj) == 1,'Can only set a single object - not a vector of objects')
            obj.y = (obj.x+eps(realmin))./(obj.x+eps(realmin)).*(obj.z+eps(realmin))./(obj.z+eps(realmin)).*y;
        end
        
        function obj = setZ(obj,z)
            %SETZ  Set the z-value of the object
            % obj = setZ(obj,z)
            assert(numel(obj) == 1,'Can only set a single object - not a vector of objects')
            obj.z = (obj.x+eps(realmin))./(obj.x+eps(realmin)).*(obj.y+eps(realmin))./(obj.y+eps(realmin)).*z;
        end
        
         %% Object operations
        function obj = getNpts(obj,I)
            %GETNPTS   Returns a subset of the points.
            % obj = getNpts(obj,I) returns a Pnt3D object only the points 
            % in the requested indexes
            %
            % Inputs
            % - obj:  Pnt3D object
            % - I:    Int array of required indexes 
            %
            % Outputs
            % - obj:  Pnt3D object
            %
            % Dependencies
            % - 
            %
            % Created: 2019-05-06, Dirk de Villiers
            % Updated: 2019-05-07, Dirk de Villiers
            %
            % Tested : Matlab R2018b, Dirk de Villiers
            %  Level : 2 
            %   File : testScript_Pnt3D.m
            %
            % Example
            %   [x,y,z] = deal(zeros(3,1));
            %   p = Pnt3D(x,y,z);
            %   p1 = p.getNpts([1;3])
            
            assert(numel(obj) == 1,'Can only get N points for a single object - not a vector of objects')
            obj = Pnt3D(obj.x(I),obj.y(I),obj.z(I));
        end
        
        function X = pointMatrix(obj)
            %POINTMATRIX   Returns a matrix of the points.
            % X = pointMatrix(obj) returns a matrix of size [3 x numel(obj.x)] 
            % containing the [x;y;z] points in its rows
            %
            % Inputs
            % - obj:  Pnt3D object
            %
            % Outputs
            % - X:    Matrix of points with rows [obj.x;obj.y;obj.z]
            %
            % Dependencies
            % - 
            %
            % Created: 2019-05-06, Dirk de Villiers
            % Updated: 2019-05-07, Dirk de Villiers
            %
            % Tested : Matlab R2018b, Dirk de Villiers
            %  Level : 2 
            %   File : testScript_Pnt3D.m
            %
            % Example
            %   [x,y,z] = deal(1:5);
            %   p = Pnt3D(x,y,z);
            %   X = p.pointMatrix
            
            assert(numel(obj) == 1,'Can only get point matrix for a single object - not a vector of objects')
            X = [obj.x(:),obj.y(:),obj.z(:)].';
        end
        
        %% Overloaded methods
        function obj = plus(obj,obj2)
            if numel(obj) == 1 && numel(obj2) == 1
                obj.x = obj.x+obj2.x;
                obj.y = obj.y+obj2.y;
                obj.z = obj.z+obj2.z;
            else
                assert(numel(obj)==numel(obj2),'Pnt3D object matrices must be the same size')
                o1 = obj.fuse;
                o2 = obj2.fuse;
                op = o1 + o2;
                obj = op.split;
            end
        end
        
        function obj = minus(obj,obj2)
            if numel(obj) == 1 && numel(obj2) == 1
                obj.x = obj.x-obj2.x;
                obj.y = obj.y-obj2.y;
                obj.z = obj.z-obj2.z;
            else
                assert(numel(obj)==numel(obj2),'Pnt3D object matrices must be the same size')
                o1 = obj.fuse;
                o2 = obj2.fuse;
                op = o1 - o2;
                obj = op.split;
            end
        end
        
        function obj = mean(obj)
            if length(obj) > 1, obj = fuse(obj); end
            obj.x = mean(obj.x);
            obj.y = mean(obj.y);
            obj.z = mean(obj.z);
        end
        
        function S = size(obj,posReturn)
            if nargin == 1, posReturn = []; end
            if isempty(posReturn)
                if numel(obj) > 1 || isempty(obj)
                    S = builtin('size',obj);
                else
                    S = size(obj.x);
                end
            else
                if numel(obj) > 1 || isempty(obj)
                    S = builtin('size',obj,posReturn);
                else
                    S = size(obj.x,posReturn);
                end
            end
        end
        
        function obj = unique(obj)
            if numel(obj) == 1 
                xyz = obj.pointMatrix;
                xyzu = unique(xyz.','rows','stable');
                obj.x = xyzu(:,1).';
                obj.y = xyzu(:,2).';
                obj.z = xyzu(:,3).';
            else
                assert(numel(obj)==numel(obj2),'Pnt3D object matrices must be the same size')
                o1 = obj.fuse;
                op = o1.unique;
                obj = op.split;
            end
        end
        
        function obj = transpose(obj)
            obj.x = obj.x.';
            obj.y = obj.y.';
            obj.z = obj.z.';
        end
%         function B = isequal(obj1,obj2)
%             tol = eps;
%             D = obj1-obj2;
%             B = all(all([abs(D.x),abs(D.y),abs(D.z)] < tol));
%             B = B && all(size(obj1) == size(obj2));
%         end

        %% Object element operations    
        function obj = scale(obj,scaleVal)
            %SCALE   Scale the components of the points
            % obj = scale(obj,scaleVal) scales the [x,y,z] components of
            % the object by scaleVal, which can be a scalar or a 3 element
            % vector [xScale,yScale,zScale].
            %
            % Inputs
            % - obj:  Pnt3D object
            % - scaleVal: Scaling factor. Can be scalar or numel(scaleVal) == 3
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
            %  Level : 2 
            %   File : testScript_Pnt3D.m
            %
            % Example
            %   [x,y,z] = deal(1:5);
            %   p = Pnt3D(x,y,z);
            %   scaleVal = [3,3,0.5].';
            %   ps = p.scale(scaleVal)
            
            if isscalar(scaleVal)
                scaleVal = ones(3,1).*scaleVal;
            else
                assert(numel(scaleVal) == 3);
            end
            obj.x = obj.x.*scaleVal(1);
            obj.y = obj.y.*scaleVal(2);
            obj.z = obj.z.*scaleVal(3);
        end
        
        function D = distanceCart(obj1,obj2)
            %DISTANCECART   Cartesian distance between points
            % D = distanceCart(obj1,obj2) returns the cartesian distance
            % between all the points in obj1 and obj2. The objects must
            % be of the same size, or either one can be scalar.
            %
            % Inputs
            % - obj1, obj2:  Pnt3D objects, must be of the same size or scalar.
            %
            % Outputs
            % - D:  Distance between all the points in the objects in m
            %
            % Dependencies
            % - 
            %
            % Created: 2019-05-07, Dirk de Villiers
            % Updated: 2019-05-07, Dirk de Villiers
            %
            % Tested : Matlab R2018b, Dirk de Villiers
            %  Level : 2 
            %   File : testScript_Pnt3D.m
            %
            % Example
            %   p = Pnt3D(1,1,1);
            %   [x,y,z] = deal([1:5]);   
            %   p1 = Pnt3D(x,y,z);
            %   pd = distanceCart(p,p1)
            
            objD = obj1 - obj2;
            D = objD.r;
        end
        
         function obj = addPoints(obj,obj2)
            %ADDPOINTS   Adds another Pnt3D object to the points.
            % obj = addPoints(obj,obj2) adds the Pnt3D object obj2 to the points.  
            % Vector must have the same number of columns as the points, 
            % or be scalar.  Same as plus, but with an explicit input
            % vector of points.
            %
            % Inputs
            % - obj, obj2:  Pnt3D objects
            %
            % Outputs
            % - obj:  Pnt3D object
            %
            % Dependencies
            % - 
            %
            % Created: 2020-03-10, Ridalise Louw
            % Updated: 2020-03-10, Ridalise Louw
            %
            % Tested : // TO DO
            %  Level :  
            %   File : testScript_Pnt3D.m
            %
            % Example
            %   [x,y,z] = deal(1:5);
            %   p = Pnt3D(x,y,z);
            %   [x2,y2,z2] = deal(6:8);
            %   p2 = Pnt3D(x2,y2,z2);
            %   p = p.addPoints(p2)
            obj = Pnt3D([obj.x obj2.x],[obj.y obj2.y],[obj.z obj2.z]);
         end
        
        function obj = addVect(obj,V)
            %ADDVECT   Adds a [3xN] vector to the points.
            % obj = addVect(obj,V) adds the [3xN] vector V to the points.  
            % Vector must have the same number of columns as the points, 
            % or be scalar.  Same as plus, but with an explicit input
            % vector of points.
            %
            % Inputs
            % - obj:  Pnt3D object
            % - V:    Vector of points to be added [3xN], where N = {1 | size(obj,2)} 
            %
            % Outputs
            % - obj:  Pnt3D object
            %
            % Dependencies
            % - 
            %
            % Created: 2019-05-06, Dirk de Villiers
            % Updated: 2019-05-07, Dirk de Villiers
            %
            % Tested : Matlab R2018b, Dirk de Villiers
            %  Level : 2 
            %   File : testScript_Pnt3D.m
            %
            % Example
            %   [x,y,z] = deal(1:5);
            %   p = Pnt3D(x,y,z);
            %   V = [3,3,3].';
            %   p1 = p.addVect(V)
            if numel(V) == 3, V = V(:); end
            if size(V,1)~=3, error('V must have 3 rows [x,y,z]'); end
            P = Pnt3D(V(1,:),V(2,:),V(3,:));
            obj = obj + P;
        end
        
        function obj = changeBase(obj,coor_new,coor_base)
            %CHANGEBASE   Changes the base coordinate system.
            % obj = changeBase(obj,coor_new,coor_base) changes the base
            % coordinate system in which the object points where defined
            % from coor_base to coor_new. If only 2 parameters are
            % provided, coor_base is assumed to be the global
            % coordinate system. This method is often used to get points in
            % the global coordinate system that were defined in some other
            % local system. 
            %
            % Inputs
            % - obj:        Pnt3D object
            % - coor_new:   New CoordinateSystem object 
            % - coor_base:  Original/Base CoordinateSystem object 
            %
            % Outputs
            % - obj:  Pnt3D object
            %
            % Dependencies
            % - 
            %
            % Created: 2019-05-06, Dirk de Villiers
            % Updated: 2019-05-07, Dirk de Villiers
            %
            % Tested : Matlab R2018b, Dirk de Villiers
            %  Level : 2 
            %   File : testScript_Pnt3D.m
            %
            % Example
            %   % Defined in a shifted and rotated coordinate system, and
            %   % plotted in the global coordinate system
            %   [x,y] = deal(1:5);
            %   z = 0;
            %   coorBase = CoordinateSystem;
            %   coorBase = coorBase.translate([0,0,5]);
            %   coorBase = coorBase.rotY(deg2rad(45));
            %   coorNew = CoordinateSystem;
            %   p = Pnt3D(x,y,z);
            %   p1 = p.changeBase(coorNew,coorBase);
            %   p.plot('marker','o'), hold on
            %   p1.plot('marker','.')
            %   coorBase.plot, coorNew.plot
            
            if nargin == 2
                coor_base = CoordinateSystem();
            end
            
            cG = CoordinateSystem;
            % Go from base to global
            if ~isequal(coor_base,cG)
                Ubase = pointMatrix(obj);
                % Rotate the points in the origin reference
                Q = dirCosine(cG,coor_base);
                Uprime = Q\Ubase;
                % Move to new coordinate base
                Uglob = Uprime + coor_base.origin.pointMatrix;
            else
                Uglob = pointMatrix(obj);
            end
            % Go from global to coorNew
            if ~isequal(coor_new,cG)
                % Move points to new coordinate origin reference
                U = Uglob - coor_new.origin.pointMatrix;
                % Rotate the points in the origin reference
                Q = dirCosine(coor_new,cG);
                Uprime = Q\U;
            else
                Uprime = Uglob;
            end
            
            % Make the object
            obj.x = reshape(Uprime(1,:),size(obj));
            obj.y = reshape(Uprime(2,:),size(obj));
            obj.z = reshape(Uprime(3,:),size(obj));
        end

        %% Plotting
        function plotVect(obj,V,varargin)
            %PLOTVECT   Plots a direction vector at the points.
            % plotVect(obj,V,varargin) plots a direction vector [3 x N]
            % at each point in the object obj on a 3-D graph. The
            % direction vector can contain only one direction, or one for
            % each of the points in the object. Similarly, the object can
            % contain only one point, or the same number as the number of
            % vectors.
            %
            % Inputs
            % - obj: Pnt3D object
            % - V:    Vector of points to be added [3xN], where N = {1 | size(obj,2)} 
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
            % Created: 2019-05-07, Dirk de Villiers
            % Updated: 2019-05-07, Dirk de Villiers
            %
            % Tested : Matlab R2018b, Dirk de Villiers
            %  Level : 2
            %   File : testScript_Pnt3D.m
            %
            % Example
            %   [x,y] = deal(1:5);
            %   z = 0; 
            %   p = Pnt3D(x,y,z);
            %   V = [0.2,2,0.7].';
            %   p.plotVect(V,'lineStyle','-','lineWidth',2), hold on
            %   p.plot('marker','*')
            %
            %   p1 = Pnt3D;
            %   V = [1:5;zeros(1,5);ones(1,5)];
            %   p1.plotVect(V,'lineStyle','-','lineWidth',2), hold on
            %   p1.plot('marker','*')
            
            parseobj = inputParser;
            parseobj.FunctionName = 'plotLines';
            
            typeValidationObj = @(x) validateattributes(x,{'Pnt3D'},{'numel',1},'plot','obj',1);
            addRequired(parseobj,'obj',typeValidationObj);

            typeValidationObj = @(x) validateattributes(x,{'double'},{'nrows',3},'plot','V',1);
            addRequired(parseobj,'V',typeValidationObj);

            typeValidationLineStyle = @(x) validateattributes(x,{'char'},{},'plot','lineStyle');
            addParameter(parseobj,'lineStyle','-',typeValidationLineStyle);
            
            typeValidationLineColor = @(x) validateattributes(x,{'char','double'},{},'plot','lineColor');
            addParameter(parseobj,'lineColor','k',typeValidationLineColor);
            
            typeValidationLineWidth = @(x) validateattributes(x,{'double'},{'real','positive'},'plot','lineWidth');
            addParameter(parseobj,'lineWidth',1,typeValidationLineWidth);
            
            parse(parseobj, obj, V, varargin{:});
            
            lineStyle = parseobj.Results.lineStyle;
            lineColor = parseobj.Results.lineColor;
            lineWidth = parseobj.Results.lineWidth;
            
            % Check the input sizes
            Nv = size(V,2);
            if ~isscalar(obj)
                x1 = obj.x(:).';
                y1 = obj.y(:).';
                z1 = obj.z(:).';
                if Nv == 1
                    V = repmat(V,1,length(x1));
                else
                    assert(max(size(obj.x)) == Nv,'There should be the same number of vectors as points')
                end
            else
                x1 = repmat(obj.x(:).',1,Nv);
                y1 = repmat(obj.y(:).',1,Nv);
                z1 = repmat(obj.z(:).',1,Nv);
            end
            
            x2 = x1 + V(1,:);
            y2 = y1 + V(2,:);
            z2 = z1 + V(3,:);
            plot3([x1;x2],[y1;y2],[z1;z2],'linestyle',lineStyle,...
                'color',lineColor,'lineWidth',lineWidth);
        end
        
        function plotLines(obj1,obj2,varargin)
            %PLOTLINES   Plots lines between two objects.
            % plotLines(obj1,obj2,varargin) plots lines between all the 
            % points in obj1 and obj2 on a 3-D graph. Both objects must be
            % the same size, or obj1 can be scalar
            % There are several options for the plotting style,
            % all a subset of the standard MATLAB plot3 options
            %
            % Inputs
            % - obj1: Pnt3D object
            % - obj2: Pnt3D object
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
            % Created: 2019-05-07, Dirk de Villiers
            % Updated: 2019-05-07, Dirk de Villiers
            %
            % Tested : Matlab R2018b, Dirk de Villiers
            %  Level : 2
            %   File : testScript_Pnt3D.m
            %
            % Example
            %   p1 = Pnt3D;
            %   x = -4:2:4;
            %   [y,z] = deal(1:5);
            %   y = y.^2 - 2;
            %   p2 = Pnt3D(x,y,z);
            %   plotLines(p1,p2,'lineStyle','--','lineWidth',1.5,'lineColor',[0.5,0.5,0.5]), hold on
            %   p1.plot('marker','*')
            %   p2.plot('marker','o')
            
            parseobj = inputParser;
            parseobj.FunctionName = 'plotLines';
            
            typeValidationObj = @(x) validateattributes(x,{'Pnt3D'},{'numel',1},'plot','obj1',1);
            addRequired(parseobj,'obj1',typeValidationObj);

            typeValidationObj = @(x) validateattributes(x,{'Pnt3D'},{'numel',1},'plot','obj2',1);
            addRequired(parseobj,'obj2',typeValidationObj);

            typeValidationLineStyle = @(x) validateattributes(x,{'char'},{},'plot','lineStyle');
            addParameter(parseobj,'lineStyle','-',typeValidationLineStyle);
            
            typeValidationLineColor = @(x) validateattributes(x,{'char','double'},{},'plot','lineColor');
            addParameter(parseobj,'lineColor','k',typeValidationLineColor);
            
            typeValidationLineWidth = @(x) validateattributes(x,{'double'},{'real','positive'},'plot','lineWidth');
            addParameter(parseobj,'lineWidth',1,typeValidationLineWidth);
            
            parse(parseobj, obj1, obj2, varargin{:});
            
            lineStyle = parseobj.Results.lineStyle;
            lineColor = parseobj.Results.lineColor;
            lineWidth = parseobj.Results.lineWidth;
            
            % Check the input sizes
            No1 = max(size(obj1.x));
            No2 = max(size(obj2.x));
            if ~isscalar(obj1) && ~isscalar(obj2)
                assert(numel(obj1.x) == numel(obj2.x),'The two points objects should be the same size')
            end
            if ~isscalar(obj1)
                x1 = obj1.x(:).';
                y1 = obj1.y(:).';
                z1 = obj1.z(:).';
            else
                x1 = repmat(obj1.x(:).',1,No2);
                y1 = repmat(obj1.y(:).',1,No2);
                z1 = repmat(obj1.z(:).',1,No2);
            end
            if ~isscalar(obj2)
                x2 = obj2.x(:).';
                y2 = obj2.y(:).';
                z2 = obj2.z(:).';
            else
                x2 = repmat(obj2.x(:).',1,No1);
                y2 = repmat(obj2.y(:).',1,No1);
                z2 = repmat(obj2.z(:).',1,No1);
            end
            plot3([x1;x2],[y1;y2],[z1;z2],'linestyle',lineStyle,...
                'color',lineColor,'lineWidth',lineWidth);
        end
        
        function plot(obj,varargin)
            %PLOT   Plots a Pnt3D object.
            % plot(obj,varargin) plots all the points in the object on a
            % 3-D graph. There are several options for the plotting style,
            % all a subset of the standard MATLAB plot3 options
            %
            % Inputs
            % - obj: Pnt3D object
            % * Arbitrary number of pairs of arguments: ...,keyword,value,... where
            %   keywords and values are from the sets
            %   -- marker:      Standard MATLAB marker chars ('.')
            %   -- markerEdgeColor: Standard MATLAB color representations ('k')
            %   -- markerFaceColor: Standard MATLAB color representations ('none')
            %   -- markerSize:  Positive real number (10)
            %   -- lineStyle:  Standard MATLAB linestyle chars ('none')
            %   -- lineColor: Standard MATLAB color representations ('k')
            %   -- lineWidth:  Positive real number (1)
            %
            % Outputs
            % - []
            %
            % Created: 2019-05-07, Dirk de Villiers
            % Updated: 2020-03-23, Dirk de Villiers
            %
            % Tested : Matlab R2018b, Dirk de Villiers
            %  Level : 2
            %   File : testScript_Pnt3D
            %
            % Example
            %   p = Pnt3D;
            %   p.plot
            %
            %   [x,y] = deal(1:5);
            %   z = 0; 
            %   p = Pnt3D(x,y,z);
            %   p.plot('marker','o','lineStyle','-','lineWidth',2)
            
            parseobj = inputParser;
            parseobj.FunctionName = 'plot';
            
            typeValidationObj = @(x) validateattributes(x,{'Pnt3D'},{},'plot','obj',1);
            addRequired(parseobj,'obj',typeValidationObj);
            
            typeValidationMarker = @(x) validateattributes(x,{'char'},{},'plot','marker');
            addParameter(parseobj,'marker','.',typeValidationMarker);
            
            typeValidationMarkerEdgeColor = @(x) validateattributes(x,{'char','double'},{},'plot','markerEdgeColor');
            addParameter(parseobj,'markerEdgeColor','k',typeValidationMarkerEdgeColor);
            
            typeValidationMarkerFaceColor = @(x) validateattributes(x,{'char','double'},{},'plot','markerFaceColor');
            addParameter(parseobj,'markerFaceColor','none',typeValidationMarkerFaceColor);
            
            typeValidationMarkerSize = @(x) validateattributes(x,{'double'},{'real','positive'},'plot','markerSize');
            addParameter(parseobj,'markerSize',10,typeValidationMarkerSize);
            
            typeValidationLineStyle = @(x) validateattributes(x,{'char'},{},'plot','lineStyle');
            addParameter(parseobj,'lineStyle','none',typeValidationLineStyle);
            
            typeValidationLineColor = @(x) validateattributes(x,{'char','double'},{},'plot','lineColor');
            addParameter(parseobj,'lineColor','k',typeValidationLineColor);
            
            typeValidationLineWidth = @(x) validateattributes(x,{'double'},{'real','positive'},'plot','lineWidth');
            addParameter(parseobj,'lineWidth',1,typeValidationLineWidth);
            
            parse(parseobj, obj, varargin{:});
            
            marker = parseobj.Results.marker;
            markerEdgeColor = parseobj.Results.markerEdgeColor;
            markerFaceColor = parseobj.Results.markerFaceColor;
            markerSize = parseobj.Results.markerSize;
            lineStyle = parseobj.Results.lineStyle;
            lineColor = parseobj.Results.lineColor;
            lineWidth = parseobj.Results.lineWidth;
            
            if numel(obj) > 1, obj = obj.fuse; end
            
            plot3(obj.x(:),obj.y(:),obj.z(:),'linestyle',lineStyle,...
                'color',lineColor,'lineWidth',lineWidth,...
                'marker',marker,'markerEdgeColor',markerEdgeColor,...
                'markerFaceColor', markerFaceColor,'markerSize',markerSize), grid on
            xlabel('x (m)')
            ylabel('y (m)')
            zlabel('z (m)')
        end
        
        %% Splitting and fusing
        function obj = split(obj1)
            %SPLIT splits the single object in a vector
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
            
            xi = reshape([obj1.x],size(obj1));
            yi = reshape([obj1.y],size(obj1));
            zi = reshape([obj1.z],size(obj1));
            assert(numel(xi) == numel(obj1),'Each object in the input array may only have a single element');
            obj = Pnt3D(xi,yi,zi);
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
        % This method is not used yet, but might be handy later
        function [objE1,objE2] = expandScalars(obj1,obj2)
            S1 = size(obj1);
            S2 = size(obj2);
            if isscalar(obj1)
                E1x = repmat(obj1.x,S2(1),S2(2));
                E1y = repmat(obj1.y,S2(1),S2(2));
                E1z = repmat(obj1.z,S2(1),S2(2));
                objE1 = Pnt3D(E1x,E1y,E1z);
            else
                objE1 = obj1;
            end
            if isscalar(obj2)
                E2x = repmat(obj2.x,S1(1),S1(2));
                E2y = repmat(obj2.y,S1(1),S1(2));
                E2z = repmat(obj2.z,S1(1),S1(2));
                objE2 = Pnt3D(E2x,E2y,E2z);
            else
                objE2 = obj2;
            end
        end
        
    end
    
end
