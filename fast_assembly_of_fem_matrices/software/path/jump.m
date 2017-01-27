function bval = jump( type, param, points, reverse )

%JUMP   A jump function
%   This function creates a boolean jump function in 2D.
%
% SYNTAX:  bval = jump( type, param, points )
%
% IN:   type      STR     'none', 'x-axis', 'y-axis', 'diagonal'
%       param     NUM     where the jump is made
%       points    2xM     the points in 2D
%       reverse   STR     'reverse'=reverse the boolean output
%
% OUT:  bval      M x 1   boolean function, where the points on the left
%                         (or above) have zero boolean values, others one
%                         (or the other way around, if the input reverse
%                         is given the string value 'reverse')
%
% EXAMPLE:
% 1)  bval = jump( 'none', ~, ~)
%     there is NO jump, so returns only ones
% 2)  bval = jump( 'x-axis', 1/4, points)
%     returns 'bval' such that the all points left of x = 1/4 have boolean
%     marking 0, and others 1.
% 3)  bval = jump( 'y-axis', -1/2, points,'reverse')
%     returns 'bval' such that the all points below of y = -1/2 have
%     boolean marking 0, and others 1.
% 4)  bval = jump( 'diagonal', ~, points)
%     returns 'bval' such that the all points above of y = x have
%     boolean marking 0, and others 1.
%


switch type
    
    case 'none'
        bval = boolean( ones(size(points,1),1) );
    
    case 'x-axis'
        bval = points(:,1) >= param ;
        
    case 'y-axis'
        bval = points(:,2) >= param ;
        
	case 'diagonal'
        bval = points(:,2) >= points(:,1) ;
        
    otherwise
        error('JUMP: The type is not recognized.')
    
end

if ( nargin > 3 && strcmp(reverse,'reverse') )
    bval = ~bval;
end