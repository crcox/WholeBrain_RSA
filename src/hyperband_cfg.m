function [n_i, r_i, n, r, B, s_max] = hyperband_cfg(max_iter, eta)
logeta = @(x) log(x)/log(eta);
s_max = floor(logeta(max_iter));
B = double((s_max+1)*max_iter);

n = zeros(s_max + 1, 1);
r = zeros(s_max + 1, 1);
n_i = cell(s_max + 1, 1);
r_i = cell(s_max + 1, 1);

x = 0;
for s = double(s_max:-1:0)
    x = x + 1;
    n(x) = ceil( ((B/max_iter)/(s+1)) * (eta.^s) );
    r(x) = floor(max_iter*eta.^(-s));

    n_i{x} = zeros(1, s+1);
    r_i{x} = zeros(1, s+1);
    y = 0; 
    for i = 0:(s)
        y = y + 1;
        n_i{x}(y) = floor(n(x)*eta.^(-i));
        r_i{x}(y) = floor(r(x)*eta.^(i));
    end
end