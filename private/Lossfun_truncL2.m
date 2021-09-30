function [l, err] = Lossfun_truncL2(x, mu)
    l = double( x < mu );
%     l = double( x/2 < mu );
%     if nargout > 1
%         err = mu * (1 - l);
%     end
%     l = double( x.^2 < mu );
%     if nargout > 1
%         err = mu * (1 - l);
%     end
end