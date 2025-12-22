function [z, x] = badpdf(x0, m)

if nargin < 2
    [z, e] = histcounts(x0);
else
    [z, e] = histcounts(x0, m);
end

x = e(1:end-1) + diff(e) / 2;
z = z / trapz(x, z);

if nargout < 1
    plot(x, z, 'k.', 'LineWidth', 1.5);
    title('Probability Density Function');
end
