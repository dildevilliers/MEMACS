classdef CoordinateSystem
    %COORDINATESYSTEM Class to describe a cartesian coordinate system
    % Cartesian coordinate systems, of arbitrary location and orientation
    % is described and handled by objects of this class. Several methods
    % for translation, rotation, and plotting is implemented herein.
    
    properties
      origin(1,1) Pnt3D = Pnt3D(0,0,0) % Coordinate system origin position
      base = [] % Coordinate system in which this one is defined. Empty implies the global system.
    end
   
    properties (SetAccess = private)
        x_axis = [1;0;0] % x-axis direction unit vector
        y_axis = [0;1;0] % y-axis direction unit vector
    end
    
    properties (Dependent = true)
        z_axis          % z-axis direction unit vector
    end
    
    methods
       function obj = CoordinateSystem(origin,x_axis,y_axis,base)
           %COORDINATESYSTEM class constructor
           % obj = CoordinateSystem(origin,x_axis,y_axis,base) Can be an
           % empty constructor which returns the standard cartesian
           % coordinate system.  The input direction vectors need not be
           % unit vectors - they are normalised internally if required.
           %
           % Inputs
           % - origin:  Pnt3D object describing the origin position
           % - x_axis:  [3x1] unit vector of the x-axis direction [ux;uy;uz]
           % - y_axis:  [3x1] unit vector of the y-axis direction [ux;uy;uz]
           % - base:    CoordinateSystem object describing the base system
           %            in which the current object is defined. Empty
           %            implies global coordinates.
           %
           % Outputs
           % - obj:  CoordinateSystem object
           %
           % Dependencies
           % -
           %
           % Created: 2019-05-09, Dirk de Villiers
           % Updated: 2019-05-09, Dirk de Villiers
           %
           % Tested : Matlab R2018b, Dirk de Villiers
           %  Level : 2
           %   File : testScript_CoordinateSystem.m
           %
           % Example
           %   C0 = CoordinateSystem;
           %   C0.plot, hold on
           %   C = CoordinateSystem(Pnt3D(1,2,3),[1;0;0],[0;-2;0],[])
           %   C.plot
           
           if nargin == 0
           elseif nargin == 1
               obj.origin = origin;
           elseif nargin == 2
               obj.origin = origin;
               obj.x_axis = x_axis;
           elseif nargin == 3
               obj.origin = origin;
               obj.x_axis = x_axis;
               obj.y_axis = y_axis;
           else
               obj.origin = origin;
               obj.x_axis = x_axis;
               obj.y_axis = y_axis;
               if isempty(base)
                   obj.base = [];
               else
                   if ~(isempty(base.base))
                       error('base coordinate system must be based on the global coordinates - use getInGlobal to get base coordinates in global coordinates')
                   else
                       obj.base = base;
                   end
               end
           end
           % Check for valid inputs
           nX = norm(obj.x_axis);
           nY = norm(obj.y_axis);
           z = cross(obj.x_axis,obj.y_axis)/(nX*nY);
           if abs(norm(z)-1) > 1e-10
               error('x_axis and y_axis must be orthogonal');
           end
           % Set the z-axis direction
           % Normalise the unit vectors
           obj = obj.normAxis;
       end
       
       %% Setters
       function z_axis = get.z_axis(obj)
            z_axis = cross(obj.x_axis,obj.y_axis);
       end
       
       function obj = set2Base(obj)
           %SET2BASE changes the current system to its base system
           % obj = set2Base(obj)
           %
           % Inputs
           % - obj:  CoordinateSystem object
           %
           % Outputs
           % - obj:  CoordinateSystem object
           %
           % Dependencies
           % -
           %
           % Created: 2019-05-09, Dirk de Villiers
           % Updated: 2019-05-09, Dirk de Villiers
           %
           % Tested : Matlab R2018b, Dirk de Villiers
           %  Level : 2
           %   File : testScript_CoordinateSystem.m
           %
           % Example
           % Cb = CoordinateSystem(Pnt3D(1,2,3),[1;0;0],[0;-2;0],[]);
           % C = CoordinateSystem(Pnt3D(1,-1,0),[1;0;0],[0;1;0],Cb);
           % C = C.set2Base;
           % Cb.plot, hold on
           % C.plot
           
           if isa(obj.base,'CoordinateSystem')
               obj = obj.base;
           end
       end
       
       %% Testers
       function  B = isequal(coor1,coor2,tol)
           %ISEQUAL tests equality between objects
           % B = isequal(coor1,coor2,tol) tests if the CoordinateSystem
           % objects defined by coor1 and coor2 are equal to within a
           % specified (optional) tolerance tol.  Ignores the bases - these
           % can be tested for equality seperately and iteratively if
           % required...
           %
           % Inputs
           % - coor1, coor2:   CoordinateSystem objects
           % - tol:   Tolerance scalar in m (1e-10);
           %
           % Outputs
           % - B:  Boolean result
           %
           % Dependencies
           % -
           %
           % Created: 2019-05-09, Dirk de Villiers
           % Updated: 2019-05-09, Dirk de Villiers
           %
           % Tested : Matlab R2018b, Dirk de Villiers
           %  Level : 2
           %   File : testScript_CoordinateSystem.m
           %
           % Example
           % C0 = CoordinateSystem;
           % C1 = CoordinateSystem(Pnt3D(1,0,0));
           % B = isequal(C0,C1,1e-12)
           % Be = isequal(C0,C0)
           
           if nargin == 2
               tol = 1e-10;
           end
           BO = isequal(coor1.origin,coor2.origin);
           Bx = all(abs(coor1.x_axis - coor2.x_axis) < tol);
           By = all(abs(coor1.y_axis - coor2.y_axis) < tol);
           B = BO && Bx && By;
       end
       
       %% Translation
       function obj = translate(obj,delta)
           %TRANSLATE translates the object by distance delta
           % obj = translate(obj,delta) shifts the coordinate system origin
           % by a distance delta = [delX,delY,delZ];
           %
           % Inputs
           % - obj:     CoordinateSystem object
           % - delta:   Vector with numel(delta)=3 by which to shift the origin in m  
           %
           % Outputs
           % - obj:  CoordinateSystem object
           %
           % Dependencies
           % -
           %
           % Created: 2019-05-09, Dirk de Villiers
           % Updated: 2019-05-09, Dirk de Villiers
           %
           % Tested : Matlab R2018b, Dirk de Villiers
           %  Level : 2
           %   File : testScript_CoordinateSystem.m
           %
           % Example
           % C0 = CoordinateSystem;
           % Cd = C0.translate([1,2,-0.5]);
           % C0.plot, hold on
           % Cd.plot
           
           obj.origin = addVect(obj.origin,delta);
       end
       
       %% Rotation
       function obj = rotX(obj,angleRadians)
           %ROTX rotates the object around the x-axis
           % obj = rotX(obj,angleRadians) rotates the object around the
           % x-axis by an angle of angleRadians in radians
           %
           % Inputs
           % - obj:         CoordinateSystem object
           % - angRadians:  Rotation angle in radians
           %
           % Outputs
           % - obj:  CoordinateSystem object
           %
           % Dependencies
           % -
           %
           % Created: 2019-05-09, Dirk de Villiers
           % Updated: 2019-05-09, Dirk de Villiers
           %
           % Tested : Matlab R2018b, Dirk de Villiers
           %  Level : 2
           %   File : testScript_CoordinateSystem.m
           %
           % Example
           % C0 = CoordinateSystem;
           % Cr = C0.rotX(deg2rad(45));
           % C0.plot, hold on
           % Cr.plot
           
           obj.x_axis = rotx(obj.x_axis,angleRadians);
           obj.y_axis = rotx(obj.y_axis,angleRadians);
       end
       
       function obj = rotY(obj,angleRadians)
           %ROTY rotates the object around the y-axis
           % obj = rotY(obj,angleRadians) rotates the object around the
           % y-axis by an angle of angleRadians in radians
           %
           % Inputs
           % - obj:         CoordinateSystem object
           % - angRadians:  Rotation angle in radians
           %
           % Outputs
           % - obj:  CoordinateSystem object
           %
           % Dependencies
           % -
           %
           % Created: 2019-05-09, Dirk de Villiers
           % Updated: 2019-05-09, Dirk de Villiers
           %
           % Tested : Matlab R2018b, Dirk de Villiers
           %  Level : 2
           %   File : testScript_CoordinateSystem.m
           %
           % Example
           % C0 = CoordinateSystem;
           % Cr = C0.rotY(deg2rad(45));
           % C0.plot, hold on
           % Cr.plot 
           
           obj.x_axis = roty(obj.x_axis,angleRadians);
           obj.y_axis = roty(obj.y_axis,angleRadians);
       end
       
       function obj = rotZ(obj,angleRadians)
           %ROTZ rotates the object around the z-axis
           % obj = rotZ(obj,angleRadians) rotates the object around the
           % z-axis by an angle of angleRadians in radians
           %
           % Inputs
           % - obj:         CoordinateSystem object
           % - angRadians:  Rotation angle in radians
           %
           % Outputs
           % - obj:  CoordinateSystem object
           %
           % Dependencies
           % -
           %
           % Created: 2019-05-09, Dirk de Villiers
           % Updated: 2019-05-09, Dirk de Villiers
           %
           % Tested : Matlab R2018b, Dirk de Villiers
           %  Level : 2
           %   File : testScript_CoordinateSystem.m
           %
           % Example
           % C0 = CoordinateSystem;
           % Cr = C0.rotZ(deg2rad(45));
           % C0.plot, hold on
           % Cr.plot 
           
           obj.x_axis = rotz(obj.x_axis,angleRadians);
           obj.y_axis = rotz(obj.y_axis,angleRadians);
       end
       
       function obj = rotGRASP(obj,angGRASP)
           %ROTGRASP rotates the object in GRASP angles
           % obj = rotGRASP(obj,angGRASP) rotates the object in the
           % GRASP angles [theta, phi, psi] in radians. They are 
           % essentially the spherical coordinate angles (theta and phi), 
           % and psi is the rotation around the local z-axis. See the 
           % GRASP_Technical_Description section 2.1
           %
           % Inputs
           % - obj:         CoordinateSystem object
           % - angGRASP:    Rotation angle in radians [theta,phi,psi]
           %
           % Outputs
           % - obj:  CoordinateSystem object
           %
           % Dependencies
           % -
           %
           % Created: 2019-05-09, Dirk de Villiers
           % Updated: 2019-05-09, Dirk de Villiers
           %
           % Tested : Matlab R2018b, Dirk de Villiers
           %  Level : 2
           %   File : testScript_CoordinateSystem.m
           %
           % Example
           % C0 = CoordinateSystem;
           % Cr = C0.rotGRASP(deg2rad([30,60,45]));
           % C0.plot, hold on
           % Cr.plot
           
           [th,ph,ps] = unpackGRASP(angGRASP);
           
           if numel(angGRASP) > 3
               error('input vector must be of length 3');
           end
           th_hat = obj.x_axis.*cos(th).*cos(ph) + obj.y_axis.*cos(th).*sin(ph) - obj.z_axis.*sin(th);
           ph_hat = -obj.x_axis.*sin(ph) + obj.y_axis.*cos(ph);
           %r_hat = obj.x_axis.*sin(th).*cos(ph) + obj.y_axis.*sin(th).*sin(ph) + obj.z_axis.*cos(th);
           obj.x_axis = th_hat.*cos(ph - ps) - ph_hat.*sin(ph - ps);
           obj.y_axis = th_hat.*sin(ph - ps) + ph_hat.*cos(ph - ps);
       end
       
       function obj = rotEuler(obj,angEuler)
           %ROTEULER rotates the object in Euler angles
           % obj = rotEuler(obj,angEuler) rotates the object in the
           % Euler angles [alpha, beta, gamma] in radians. They are 
           % essentially the GRASP angles with some offsets. See the 
           % GRASP_Technical_Description section 2.1, and GRASP2Euler.m
           %
           % Inputs
           % - obj:         CoordinateSystem object
           % - angEuler:    Rotation angle in radians [alpha,beta,gamma]
           %
           % Outputs
           % - obj:  CoordinateSystem object
           %
           % Dependencies
           % -
           %
           % Created: 2019-05-09, Dirk de Villiers
           % Updated: 2019-05-09, Dirk de Villiers
           %
           % Tested : Matlab R2018b, Dirk de Villiers
           %  Level : 2
           %   File : testScript_CoordinateSystem.m
           %
           % Example
           % C0 = CoordinateSystem;
           % Cr = C0.rotEuler(deg2rad([30,60,45]));
           % C0.plot, hold on
           % Cr.plot
           
           if numel(angEuler) > 3
               error('input vector must be of length 3');
           end
           obj = rotGRASP(obj,Euler2GRASP(angEuler));
       end
       
       
       %% Change of basis
       function Q = dirCosine(coor_new,coor_base)
           %DIRCOSINE calculates the direction cosine angles between coordinate systems
           % Q = dirCosine(coor_new,coor_base)) Calculates the [3x3]
           % direction cosine matrix between coordinate systems.  
           % For only one argument, the global coordinate system
           % is assumed as the second (coor_base) system. 
           % coor_new indicates the rotated system and coor_base the base system.
           % So, a point in the rotated system, which was specified in the
           % base system, is calculated as A_rotate = inv(Q)*A_base.
           % See the note: http://homepages.engineering.auckland.ac.nz/~pkel015/SolidMechanicsBooks/Part_III/Chapter_1_Vectors_Tensors/Vectors_Tensors_05_Coordinate_Transformation_Vectors.pdf
           %
           % Inputs
           % - coor_new:    CoordinateSystem object that was roted in the
           %                base coordinate system, coor_base
           % - coor_base:   Base CoordinateSystem object (global coor syst)
           %
           % Outputs
           % - Q:  [3x3] direction cosine matrix between the systems in rad
           %
           %
           % Dependencies
           % -
           %
           % Created: 2019-05-09, Dirk de Villiers
           % Updated: 2019-05-09, Dirk de Villiers
           %
           % Tested : Matlab R2018b, Dirk de Villiers
           %  Level : 2
           %   File : testScript_CoordinateSystem.m
           %
           % Example
           % C0 = CoordinateSystem;
           % Cr = C0.rotEuler(deg2rad([30,60,45]));
           % Q = dirCosine(Cr,C0)
           % [Cr.x_axis,Cr.y_axis,Cr.z_axis]
           
           if nargin == 1
               coor_base = CoordinateSystem;
           end
           Q = [dot(coor_base.x_axis,coor_new.x_axis), dot(coor_base.x_axis,coor_new.y_axis), dot(coor_base.x_axis,coor_new.z_axis);...
                dot(coor_base.y_axis,coor_new.x_axis), dot(coor_base.y_axis,coor_new.y_axis), dot(coor_base.y_axis,coor_new.z_axis);...
                dot(coor_base.z_axis,coor_new.x_axis), dot(coor_base.z_axis,coor_new.y_axis), dot(coor_base.z_axis,coor_new.z_axis)];
       end
       
       function angGRASP = getGRASPangBetweenCoors(coor1,coor0)
           %GETGRASPANGBETWEENCOORS calculates the GRASP angles between coordinate systems
           % Returns the GRASP angles (in rad) required to rotate from
           % coor0 to coor1.  That is coor1 = coor0.rotGRASP([th,ph,ps])
           % See the GRASP technical description for details on the
           % definitions.( coor0 corresponds to xyz; coor1 corresponds to
           % x1y1z1)
           % For only one argument, the global system is assumed for coor0
           %
           % Inputs
           % - coor0:   Base CoordinateSystem object
           % - coor1:   Rotated CoordinateSystem object
           %
           % Outputs
           % - angGRASP:    [th,ph,ps] GRASP angles in radians to rotate from
           %                coor0 to coor1
           %
           % Dependencies
           % -
           %
           % Created: 2019-05-09, Dirk de Villiers
           % Updated: 2019-05-09, Dirk de Villiers
           %
           % Tested : Matlab R2018b, Dirk de Villiers
           %  Level : 2
           %   File : testScript_CoordinateSystem.m
           %
           % Example
           % C0 = CoordinateSystem;
           % Cr = C0.rotGRASP(deg2rad([30,60,45]));
           % angGRASP = rad2deg(getGRASPangBetweenCoors(Cr,C0))
           
           if nargin == 1
               coor0 = CoordinateSystem();
           end
           
           x = coor0.x_axis;
           z = coor0.z_axis;
           x1 = coor1.x_axis;
           z1 = coor1.z_axis;
           
           % Get th - angles between z and z1
           th = angBetweenVectors(z,z1);
           
           % Now get ph
           % Get the vector which is x rotated by ph in the x-y plane 
           if ~isequal(abs(z),abs(z1))
               Nz = cross(z,z1);
               Nz = Nz./norm(Nz);
               x_ph = cross(Nz,z);
               x_ph = x_ph./norm(x_ph);
           else
               x_ph = x1;
           end
           ph = angBetweenVectors(x,x_ph);
           % Sort out the sign of ph - compare to the z-axis direction
           phSign = sign(dot(cross(x,x_ph),z));
           if phSign ~= 0
               ph = ph*phSign;
           end
           
           % Rotate the system to get to the intermediate coordinate system
           % xyz_prime
           coorPrime = coor0.rotGRASP([th,ph,0]);
           xp = coorPrime.x_axis;
           % Calculate psi
           ps = angBetweenVectors(xp,x1);
           % Sort out the sign of ps - compare to the z1-axis direction
           psSign = sign(dot(cross(xp,x1),z1));
           if psSign ~= 0
               ps = ps*psSign;
           end
           angGRASP = [th,ph,ps];
       end
       
       function angEuler = getEulerangBetweenCoors(coor1,coor0)
           %GETEULERANGBETWEENCOORS calculates the Euler angles between coordinate systems
           % Returns the Euler angles (in rad) required to rotate from
           % coor0 to coor1.  That is coor1 = coor0.rotEuler([alpha,beta,gamma])
           % For only one argument, the global system is assumed for coor0
           %
           % Inputs
           % - coor0:   Base CoordinateSystem object
           % - coor1:   Rotated CoordinateSystem object
           %
           % Outputs
           % - angEuler:    [alpha,beta,gamma] Euler angles in radians to 
           %                rotate from coor0 to coor1
           %
           % Dependencies
           % -
           %
           % Created: 2019-05-09, Dirk de Villiers
           % Updated: 2019-05-09, Dirk de Villiers
           %
           % Tested : Matlab R2018b, Dirk de Villiers
           %  Level : 2
           %   File : testScript_CoordinateSystem.m
           %
           % Example
           % C0 = CoordinateSystem;
           % Cr = C0.rotEuler(deg2rad([30,60,45]));
           % angEulerP = rad2deg(getEulerangBetweenCoors(Cr,C0))
           
           angEuler = GRASP2Euler(getGRASPangBetweenCoors(coor1,coor0));
       end
       
       function coorInGlobal = getInGlobal(obj)
           %GETINGLOBAL returns the object in the global coordinate system
           % coorInGlobal = getInGlobal(obj) Returns the coordinate system,
           % defined in obj.base, in the global coordinate system
           % Assume obj.base is always defined in the global coordinate system
           %
           % Inputs
           % - obj:  CoordinateSystem object
           %
           % Outputs
           % - coorInGlobal:  CoordinateSystem object with global base
           %
           % Dependencies
           % -
           %
           % Created: 2019-05-09, Dirk de Villiers
           % Updated: 2019-05-09, Dirk de Villiers
           %
           % Tested : Matlab R2018b, Dirk de Villiers
           %  Level : 2
           %   File : testScript_CoordinateSystem.m
           %
           % Example
           % C0 = CoordinateSystem;
           % C1 = C0.rotGRASP(deg2rad([20,-35,139]));
           % C0.base = C1;
           % C2 = C0.getInGlobal; % Should be the same as C1
           % C1.plot, hold on; C2.plot 

           if ~isempty(obj.base)
               oldBase = obj.base;
               newBase = CoordinateSystem();
               % Since the coordinate system is always defined in the global
               % coordinate system, with a given base (which corresponds to the
               % local global coordinate system), get the angle between the
               % new base coordinate system and the object
               graspAng = getGRASPangBetweenCoors(obj);
               % Now rotate the base by this amount
               baseRotated = oldBase.rotGRASP(graspAng);
               % And shift the origin
               Ups = changeBase(obj.origin,newBase,oldBase);
               coorInGlobal = CoordinateSystem(Ups,baseRotated.x_axis,baseRotated.y_axis);
               coorInGlobal.base = [];
           else
               coorInGlobal = obj;
           end
       end
       
       function coorOut = redefineToOtherBase(obj,newBase)
           %REDEFINETOOTHERBASE returns the object defined in a new base
           % coorOut = redefineToOtherBase(obj,newBase)) Redefines the 
           % coordinate system in obj to be relative to the new base 
           % coordinate system newBase
           %
           % Inputs
           % - obj:     CoordinateSystem object in some unwanted base
           % - newBase: CoordinateSystem object of the new base
           %
           % Outputs
           % - coorOut:  CoordinateSystem object in the new base
           %
           % Dependencies
           % -
           %
           % Created: 2019-05-09, Dirk de Villiers
           % Updated: 2019-05-09, Dirk de Villiers
           %
           % Tested : Matlab R2018b, Dirk de Villiers
           %  Level : 2
           %   File : testScript_CoordinateSystem.m
           %
           % Example
           % C0 = CoordinateSystem;
           % C1 = C0.rotGRASP(deg2rad([20,-35,139]));
           % C0.base = C1;
           % % C2, in global base, should be the same as C0, in C1 base
           % C2 = C0.redefineToOtherBase(CoordinateSystem);
           % C0.plot, hold on; C2.plot 
           
           % First get both coordinate systems in the global base
           coorIn = obj.getInGlobal;
           baseIn = newBase.getInGlobal;
           % Calculate the translation
           diffGlobal = coorIn.origin.pointMatrix - baseIn.origin.pointMatrix;
           Ox = dot(diffGlobal,baseIn.x_axis);
           Oy = dot(diffGlobal,baseIn.y_axis);
           Oz = dot(diffGlobal,baseIn.z_axis);
           % Calculate the rotation
           graspAng = getGRASPangBetweenCoors(coorIn,baseIn);
           % Build the coordinate system
           coorOut = CoordinateSystem(Pnt3D(Ox,Oy,Oz));
           coorOut = coorOut.rotGRASP(graspAng);
           coorOut.base = newBase;
       end
      
       %% Plotting
       function plot(obj,scale)
           %PLOT plots the object in global coordinates
           % plot(obj,scale) Plots the object in global coordinates, and
           % scales all three axis by a factor scale
           %
           % Inputs
           % - obj:     CoordinateSystem object 
           % - scale:   Scaling factor (real scalar) 
           %
           % Outputs
           % - 
           %
           % Dependencies
           % -
           %
           % Created: 2019-05-09, Dirk de Villiers
           % Updated: 2019-05-09, Dirk de Villiers
           %
           % Tested : Matlab R2018b, Dirk de Villiers
           %  Level : 2
           %   File : testScript_CoordinateSystem.m
           %
           % Example
           % C0 = CoordinateSystem;
           % C1 = C0.rotGRASP(deg2rad([20,-35,139]));
           % C0.base = C1;
           % C0.plot(2)
           
           if nargin  == 1
               scale = 1;
           end
           coorGlobalBase = obj.getInGlobal;
           plotLocal(coorGlobalBase,scale);
       end
       
       function plotLocal(obj,scale)
           %PLOTLOCAL plots the object in global coordinates
           % plotLocal(obj,scale) Plots the object in the local base
           % coordinates, and scales all three axis by a factor scale
           %
           % Inputs
           % - obj:     CoordinateSystem object 
           % - scale:   Scaling factor (real scalar) 
           %
           % Outputs
           % - 
           %
           % Dependencies
           % -
           %
           % Created: 2019-05-09, Dirk de Villiers
           % Updated: 2019-05-09, Dirk de Villiers
           %
           % Tested : Matlab R2018b, Dirk de Villiers
           %  Level : 2
           %   File : testScript_CoordinateSystem.m
           %
           % Example
           % C0 = CoordinateSystem;
           % C1 = C0.rotGRASP(deg2rad([20,-35,139]));
           % C0.base = C1;
           % C0.plotLocal(2)
           
           if nargin  == 1
               scale = 1;
           end
           lineWidth = 1;
           textSize = 12;
           x = obj.x_axis.*scale;
           y = obj.y_axis.*scale;
           z = obj.z_axis.*scale;
           O = [obj.origin.x;obj.origin.y;obj.origin.z];
           xAx = [O, O+x];
           yAx = [O, O+y];
           zAx = [O, O+z];
           plot3(xAx(1,:),xAx(2,:),xAx(3,:),'r','LineWidth',lineWidth), hold on
           plot3(yAx(1,:),yAx(2,:),yAx(3,:),'g','LineWidth',lineWidth)
           plot3(zAx(1,:),zAx(2,:),zAx(3,:),'b','LineWidth',lineWidth)
           text(xAx(1,2),xAx(2,2),xAx(3,2),'x','FontSize',textSize,'color','r')
           text(yAx(1,2),yAx(2,2),yAx(3,2),'y','FontSize',textSize,'color','g')
           text(zAx(1,2),zAx(2,2),zAx(3,2),'z','FontSize',textSize,'color','b')
           grid on
           axis equal
           xlabel('x-axis (m)')
           ylabel('y-axis (m)')
           zlabel('z-axis (m)')
           view([140,40])
       end
       
       %% Deprecated
       function angGRASP = getGRASPangles(coor)
           %GETGRASPANGLES deprecated: use getGRASPangBetweenCoors
           
           warning('Deprecated function: Use getGRASPangBetweenCoors instead')
           angGRASP = getGRASPangBetweenCoors(coor);
       end
       
       function angEuler = getEulerAngles(coor)
           %GETEULERANGLES deprecated: use getEulerangBetweenCoors
           
           warning('Deprecated function: Use getEulerangBetweenCoors instead')
           angEuler = GRASP2Euler(getGRASPangBetweenCoors(coor));
       end
   end
   
   methods (Access = private)
       function obj = normAxis(obj)
           % Make sure we have unit vectors
           obj.x_axis = obj.x_axis/norm(obj.x_axis);
           obj.y_axis = obj.y_axis/norm(obj.y_axis);
       end
   end
end
