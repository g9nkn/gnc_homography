function l = Lossfun_L1(x, mu)
    l = sqrt(mu) ./ (2 * abs(x));
end