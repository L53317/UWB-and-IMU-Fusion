function c = cros(a, b)
% The cross product of 3-element vectors a and b, i.e. c = a x b. It is
% a simple version and then much faster than Matlab lib-function 'cross'.
%
% Prototype: c = cros(a, b)
% Inputs: a, b - 3-element vectors
% Output: c - c = a x b
%
% See also  askew, iaskew, rv2m, rv2q, cros3.

% Copyright(c) 2009-2014, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 04/05/2014
    c = a;
    c(1) = a(2)*b(3)-a(3)*b(2);
    c(2) = a(3)*b(1)-a(1)*b(3);
    c(3) = a(1)*b(2)-a(2)*b(1);
    
%     if size(a,1)==1  % row vector
%         c = [ a(2)*b(3)-a(3)*b(2),...
%               a(3)*b(1)-a(1)*b(3),...
%               a(1)*b(2)-a(2)*b(1) ];
%     else
%         c = [ a(2)*b(3)-a(3)*b(2)
%               a(3)*b(1)-a(1)*b(3)
%               a(1)*b(2)-a(2)*b(1) ];
%     end
