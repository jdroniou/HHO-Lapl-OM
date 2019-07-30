% ----------------------------------------------------------------------
%                   Basis function on a 2D cell edge
% ----------------------------------------------------------------------
%
% Author: Daniel Anderson
%
% ----------------------------------------------------------------------
classdef EdgeBasisFunction
    %EDGEBASISFUNCTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        center              % The midpoint of the edge
        vertex              % A vertex of the edge
        hF                  % The length of the edge
        exponent            % The exponent of the basis function
    end
    
    methods
        function obj = EdgeBasisFunction(center,vertex,hF,exponent)
            %EDGEBASISFUNCTION Construct an edge basis function
            obj.center = center;
            obj.vertex = vertex;
            obj.hF = hF;
            obj.exponent = exponent;
        end
        
        function f = eval(obj,x,y)
            %eval Evaluate the basis function at the point (x,y)
            f = (dot([x,y] - obj.center, obj.vertex - obj.center) / obj.hF^2)^obj.exponent;
        end
    end
end

