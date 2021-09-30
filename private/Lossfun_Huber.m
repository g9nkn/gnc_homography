function [l, err] = Lossfun_Huber(x, mu)
    ind = x < mu;
    l = double( ind );
    l(~ind) = mu/sqrt(x(~ind));
    
    if nargout > 1
        err = x;
        err(~ind) = sqrt(mu*x(~ind)) - mu/2;
    end
end