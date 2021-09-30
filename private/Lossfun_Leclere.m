function l = Lossfun_Leclere(x, mu, sigma)
    l = ( mu ./ (mu + x.^2) ).^2;
end