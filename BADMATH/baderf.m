function v = baderf(varargin)
%BADERF   @(x, n, T) 4/T * (x*sqrt(2)/2 + sum(exp(-s^2 * ((1:n).').^2/2 * (2*pi/T)^2) .* sin((1:n).' .* x * sqrt(2) * (2*pi/T)) ./ ((1:n).' * (2*pi/T))));
%
%   v = BADERF() without input arguments returns a default error function 
%   approximation, v.
%
%   BADERF(...) without output arguments returns nothing, but plots the approximation.
%
%   v = BADERF(n) generates an n-level approximation with default inputs, returning
%   the approximation, v.
%
%   v = BADERF(x, n) generates an n-level approximation with given inputs, x, 
%   returning the approximation, v.
%
%   v = BADERF(x, n, T) generates an n-level approximation with given inputs, x,
%   and period, T, returning the approximation, v.
%
%   See also: BADCOS, BADSIN, BADGAUSS, BADCDF 
switch nargin
    case 0
        % v = baderf
        n = 10;
        x = -2*pi:0.01:2*pi;
        T = 10;
        
    case 1
        % v = baderf(n)
        if length(varargin{1}) > 1
            error('baderf:InvalidInput','With one input, it must be the number of harmonics to use, n.');
        end
        
        n = varargin{1};
        x = -2*pi:0.01:2*pi;
        T = 10;
        
    case 2
        % v = baderf(x, n)
        x = varargin{1};
        
        if length(varargin{2}) > 1
            error('baderf:InvalidInput','Input variable n is the number of harmonics and must be a single integer.');
        end
        n = varargin{2};
        T = 10;
                
    case 3
        % v = baderf(x, n, T)
        x = varargin{1};
        
        if length(varargin{2}) > 1 || length(varargin{3}) > 1
            error('baderf:InvalidInput','Other than x, input variables must be single numbers.');
        end
        n = varargin{2};
        T = varargin{3};
end

x = x * sqrt(2);
w = 2 * pi / T;
v = exp(-(w).^2 / 2) .* sin(w.*x) ./ (w);
for i = 2:n
    v = v + exp(-(w * i).^2 / 2) .* sin(w * i .* x) ./ (w * i);
end
v = (v + x/2) * 4 / T;

switch nargout
    case 0
        fig = figure("OuterPosition", [200, 200, 800, 500]);
        ax = gca(fig);
        plot(ax, x / sqrt(2), v, 'LineWidth', 1.5);
        grid on;
        axis([x(1), x(end), -1.1, 1.1]);
        legend('"erf"');
        title("Error Function Approximation");
        xlabel("input parameter (a.u.)");

        v = [];
end