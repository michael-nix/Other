function x = badrand(pdf, xrange, yrange, n)
    
if nargin < 1
    pdf = @(x) exp(-x.^2/2)/sqrt(2*pi);
    xrange = [-5, 5];
    yrange = [0, 1/sqrt(2*pi)];
    n = 2^17;
end

r = abs(diff(xrange) * diff(yrange));
x0 = xrange(1) + diff(xrange) * rand(1, floor(r*n));
y0 = yrange(1) + diff(yrange) * rand(1, floor(r*n));

idx = y0 <= pdf(x0);
x = x0(idx);

while length(x) < n
    x0 = xrange(1) + diff(xrange) * rand(1, n);
    y0 = yrange(1) + diff(yrange) * rand(1, n);
    
    idx = y0 <= pdf(x0);
    x = [x, x0(idx)]; %#ok<AGROW>
end
x = x(1:n);
