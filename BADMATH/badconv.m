function c = badconv(a, b)
%   c = BADCONV(a, b) computes the convolution of a and b, poorly.

if isrow(a)
    a = a.';
end

if isrow(b)
    b = b.';
end

M = zeros(length(b) + length(a) - 1, length(a));
for i = 1:length(a)
    M(i:(i+length(b)-1), i) = b;
end
c = M * a(end:-1:1);
