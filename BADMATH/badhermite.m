function h = badhermite(n)
%   h = BADHERMITE(n) returns the nth Hermite polynomial.

if n == 0
    h = 1;
    return;
end

h = zeros(1, n + 1);
h(end-1:end) = [-1, 0];
for i = 2:n
    idx = n+2-i:n+1;
    h(end-i:end) = [0, 0, polyder(h(idx))] + conv([-1, 0], h(idx));
end

h = (-1)^mod(n, 2) * h;