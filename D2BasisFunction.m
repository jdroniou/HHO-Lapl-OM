% ----------------------------------------------------------------------
%             2D      Basis function on a 2D cell
% ----------------------------------------------------------------------
%
% Author: Daniel Anderson
%
% ----------------------------------------------------------------------
classdef D2BasisFunction %basis function for a polynomial in dimension 2
    properties
        center              % The midpoint of the cell
        hT                  % The size of the cell (usually area/perimeter)
        c1                   % determines if we include the x component of the vector in the basis
        x1_pow               % The exponent of the x term, x component
        y1_pow               % The exponent of the y term, x component
        c2                   % determines if we include the y component of the vector in the basis
        x2_pow               % The exponent of the x term, y component
        y2_pow               % The exponent of the y term, y component
    end
    
    methods
        function obj = D2BasisFunction(center,hT,c1,x1_pow,y1_pow,c2,x2_pow,y2_pow)
            %RTBasisFunction Construct a 2D basis function
            obj.center = center;
            obj.hT = hT;
            obj.c1=c1;
            obj.x1_pow = x1_pow;
            obj.y1_pow = y1_pow;
            obj.c2=c2;
            obj.x2_pow = x2_pow;
            obj.y2_pow = y2_pow;
        end
        
        function f = eval(obj,x,y)
            %eval Evaluate the basis function at the point (x,y)
            % basis functions of the form x^ny^m
            %% Note: for vectorialised evaluation of RT functions, this works better than the latter function, 
            %% However, when solving the local systems, yield matrices with very high condition numbers (not good)
%              f = [obj.c1 .* (x) .^ obj.x1_pow .* (y) .^ obj.y1_pow; ...
%                 obj.c2 .* (x) .^ obj.x2_pow .* (y) .^ obj.y2_pow];
            % basis functions of the form (x-x_K)^n(y-y_K)^m/h^{m+n}
            %% Local systems are better defined with this basis, but can't have vectorialised implementation
            f = [obj.c1*((x - obj.center(1))/obj.hT)^obj.x1_pow * ((y - obj.center(2))/obj.hT)^obj.y1_pow; ...
                obj.c2*((x - obj.center(1))/obj.hT)^obj.x2_pow * ((y - obj.center(2))/obj.hT)^obj.y2_pow];
        end
    end
end

