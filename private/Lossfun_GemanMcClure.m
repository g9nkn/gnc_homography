function [l, err] = Lossfun_GemanMcClure(x, mu)
    l = ( mu ./ (mu + x) ).^2;
%     l = ( (2*mu)^2 ./ (2*mu + x) ).^2;
    
    if nargout > 1
        err = mu * x./(1 + x);
    end
%     l = ( mu ./ (mu + x.^2) ).^2;
%     
%     if nargout > 1
%         err = mu * (x.^2)./(1 + x.^2);
%     end
end