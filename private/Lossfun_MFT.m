function [l, err] = Lossfun_MFT(x, mu)
    l = 1 ./ ( exp(x/mu - 1) + 1 );
    if nargout > 1
        err = mu* ( 1-x+eps + (1-x+eps).*log(1-x+eps) + (x+eps).*log(x+eps) );
    end
%     z = x.^2;
%     l = 1 ./ ( exp(z/mu - 1) + 1 );
%     if nargout > 1
%         err = mu* ( 1-z+eps + (1-z+eps).*log(1-z+eps) + (z+eps).*log(z+eps) );
%     end
end