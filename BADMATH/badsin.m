function [u, x] = badsin(varargin)
%BADSIN   @(x, n) sum(1i.^((0:2:n-1)')./factorial((1:2:n)').*x.^((1:2:n)'));
%
%   [u, x] = BADSIN() without input arguments returns a default sine approximation,
%   u, along with the inputs that were used to generate it.
%
%   BADSIN(...) without output arguments returns nothing, but plots the approximation.
%
%   [u, x] = BADSIN(n) generates an n-level approximation with default inputs, returning
%   the approximation, u, and the inputs used to generate it, x.
%
%   [u, x] = BADSIN(x, n) generates an n-level approximation with given inputs, x,
%   returning the approximation, u, and the inputs used to generate it, x.
%
%   See also: BADCOS, BADGAUSS, BADCDF, BADERF
if nargin > 1 && length(varargin) < 3
    x = varargin{1};
    n = varargin{2};
elseif nargin == 1
    if length(varargin{1}) > 1
        error('If single argument, input polynomial order.');
    end
    x = linspace(-pi, pi, 100);
    n = varargin{1};
elseif nargin == 0
    x = linspace(-pi, pi, 100);
    n = 5;
end

n = (1:2:n).';
u = sum(1i.^(n-1) ./ factorial(n) .* x.^n, 1);

switch nargout
    case 0
        fig = figure("OuterPosition", [200, 200, 800, 500]);
        ax = gca(fig);
        plot(ax, x, u, 'LineWidth', 1.5);
        grid on;
        axis([x(1), x(end), -1.1, 1.1]);
        legend('"sine"');
        title("Polynomial Sine Function Approximation");
        xlabel("input parameter (a.u.)");
        ylabel("sine")

        u = [];
end