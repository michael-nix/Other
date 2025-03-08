function v = badcdf(varargin)
% @(x, n, s, T) 2/T*(x/2 + sum(exp(-s^2*((1:n).').^2/2*(2*pi/T)^2).*sin((1:n).'.*x*(2*pi/T)) ./ ((1:n).'*(2*pi/T)))) + 1/2;

switch nargin
    case 0
        % [u, x] = badgauss
        n = 10;
        x = -2*pi:0.01:2*pi;
        s = 0.5;
        T = 2 * 5 * s;
        
    case 1
        % [u, x] = badgauss(n)
        if length(varargin{1}) > 1
            error('badgauss:InvalidInput','With one input, it must be the number of harmonics to use, n.');
        end
        
        n = varargin{1};
        x = -2*pi:0.01:2*pi;
        s = 0.5;
        T = 2 * 5 * s;
        
    case 2
        % [u, x] = badgauss(x, n)
        x = varargin{1};
        
        if length(varargin{2}) > 1
            error('badgauss:InvalidInput','Input variable n is the number of harmonics and must be a single integer.');
        end
        n = varargin{2};
        s = 0.5;
        T = 2 * 5 * s;
        
    case 3
        % [u, x] = badgauss(x, n, s)
        x = varargin{1};
        
        if length(varargin{2}) > 1 || length(varargin{3}) > 1
            error('badgauss:InvalidInput','Other than x, input variables must be single numbers.');
        end
        n = varargin{2};
        s = varargin{3};
        T = 2 * 5 * s;
        
    case 4
        % [u, x] = badgauss(x, n, s, T)
        x = varargin{1};
        
        if length(varargin{2}) > 1 || length(varargin{3}) > 1 || length(varargin{4}) > 1 
            error('badgauss:InvalidInput','Other than x, input variables must be single numbers.');
        end
        n = varargin{2};
        s = varargin{3};
        T = varargin{4};
end

w = 2 * pi / T;
v = exp(-(s * w).^2 / 2) .* sin(w.*x) ./ (w);
for i = 2:n
    v = v + exp(-(s * w * i).^2 / 2) .* sin(w * i .* x) ./ (w * i);
end
v = (v + x/2) * 2 / T + 1/2;

% slower:
% n = (1:n)';
% w = 2*pi*n/T;
% v = 2/T*(x/2 + sum(exp(-(s*w).^2/2).*sin(w.*x)./w)) + 1/2;

switch nargout
    case 0
        fig = figure("OuterPosition", [200, 200, 800, 500]);
        ax = gca(fig);
        plot(ax, x, v, 'LineWidth', 1.5);
        grid on;
        axis([x(1), x(end), -0.1, 1.1]);
        legend('"CDF"');
        title("Gaussian Cumulative Density Function Approximation");
        xlabel("input parameter (a.u.)");
        ylabel("probability");

        v = [];
end