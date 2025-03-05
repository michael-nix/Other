function [u, x] = badcos(varargin)
% @(x,n) sum(1i.^((0:2:n).')./factorial((0:2:n).').*x.^((0:2:n).'),1);

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

n = (0:2:n).';
u = sum(1i.^n ./ factorial(n) .* x.^n, 1);

switch nargout
    case 0
        fig = figure("OuterPosition", [200, 200, 800, 600]);
        ax = gca(fig);
        plot(ax, x, u, 'LineWidth', 1.5);
        grid on;
        axis([min(x), max(x), -1.1, 1.1]);
        legend('"cosine"');
        title("Polynomial Cosine Function Approximation");
        xlabel("input parameter (a.u.)");
        ylabel("cosine")

        u = [];
end