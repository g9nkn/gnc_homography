function [l, err] = Lossfun_Tukey(x, mu)
    
    ind    = x < mu;
    l      = zeros(length(x),1);
    l(ind) = ( 1 - (x(ind)/mu).^2 ).^2;
    
end