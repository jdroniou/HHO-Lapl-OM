% ----------------------------------------------------------------------
%                   Basis function on a 2D cell
% ----------------------------------------------------------------------
%
% Author: Daniel Anderson
%
% ----------------------------------------------------------------------
classdef CellBasisFunction
    properties
        center              % The midpoint of the cell
        hT                  % The size of the cell (usually area/perimeter)
        x_pow               % The exponent of the x term
        y_pow               % The exponent of the y term
    end
    
    methods
        function obj = CellBasisFunction(center,hT,x_pow,y_pow)
            %CELLBASISFUNCTION Construct a cell basis function
            obj.center = center;
            obj.hT = hT;
            obj.x_pow = x_pow;
            obj.y_pow = y_pow;
        end
        
        function f = eval(obj,x,y)
            %eval Evaluate the basis function at the point (x,y)
            f = ((x - obj.center(1))/obj.hT)^obj.x_pow * ((y - obj.center(2))/obj.hT)^obj.y_pow;
        end
        
        function g = eval_grad(obj,x,y)
            %eval_grad Evaluate the gradient of the basis function at the point (x,y)
            g = [0 0];
            if obj.x_pow ~= 0
                g(1) = obj.x_pow / obj.hT * ((x - obj.center(1))/obj.hT)^(obj.x_pow-1) * ((y - obj.center(2))/obj.hT)^obj.y_pow;
            end
            if obj.y_pow ~= 0
                g(2) = ((x - obj.center(1))/obj.hT)^obj.x_pow * obj.y_pow / obj.hT * ((y - obj.center(2))/obj.hT)^(obj.y_pow-1);
            end
        end
    end
end

