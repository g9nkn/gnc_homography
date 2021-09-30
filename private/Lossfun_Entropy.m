function [l, err] = Lossfun_Entropy(x, mu)
    l = 2./(1 + exp(x/mu));
    if nargout > 1
        x2 = sqrt(x);
        err = mu * ( (x2+eps).*log(x2+eps) + (1-x2+eps).*log(1-x2+eps) );
    end
%     l = 2./(1 + exp(x.^2/mu));
%     if nargout > 1
%         err = mu * ( (x+eps).*log(x+eps) + (1-x+eps).*log(1-x+eps) );
%     end
end