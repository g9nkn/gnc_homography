function [l, err] = Lossfun_L2log(x, mu)
    l = min( mu./x, 1);
    
    if nargout > 1
        ind = (x < mu);
        err = x.*ind + (mu*(log(x/mu) + 1)).*(~ind);
    end
    
% %     l = mu./(x.^2);
%     l = min( mu./(x.^2), 1);
%     
%     if nargout > 1
%         x2  = x.^2;
%         ind = (x.^2 < mu);
%         err = (x2).*ind + (mu*(log(x2/mu) + 1)).*(~ind);
%     end
    
end