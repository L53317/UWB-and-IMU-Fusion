function c = cros(a, b)
% The cross product of 3-element vectors a and b, i.e. c = a x b. It is
% a simple version and then much faster than Matlab lib-function 'cross'.
%
% Prototype: c = cros(a, b)
% Inputs: a, b - 3-element vectors
% Output: c - c = a x b

    c = a;
    c(1) = a(2)*b(3)-a(3)*b(2); % 旨在直接快速计算，已不需要，可参考cross
    c(2) = a(3)*b(1)-a(1)*b(3);
    c(3) = a(1)*b(2)-a(2)*b(1);
end