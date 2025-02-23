function v = badcdf(x, n, s, T)
    w = 2 * pi / T;
    v = exp(-(s * w).^2 / 2) .* sin(w.*x) ./ (w);
    for i = 2:n
        v = v + exp(-(s * w * i).^2 / 2) .* sin(w * i .* x) ./ (w * i);
    end
    v = (v + x/2) * w / pi + 1/2;

    % slower:
    % n = (1:n)';
    % w = 2*pi*n/T;
    % v = 2/T*(x/2 + sum(exp(-(s*w).^2/2).*sin(w.*x)./w, 1)) + 1/2;