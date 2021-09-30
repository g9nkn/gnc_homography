function y = H_reproj_lsqnonlin(x, m1, m2, weight, options)
    
    n = size(m1,2);
    if nargin < 4
        weight = ones(1, n);
    end
    if isrow(weight)
        weight = weight';
    end
    
%     options.Algorithm = 'trust-region-reflective';

    x0 = [x(1:8);
          m1(1,:)';
          m1(2,:)'];    
    x = lsqnonlin(@(x)objfun(x,m1(1:2,:),m2(1:2,:),weight), x0, [], [], options);
    
    
    % calc optimal x2 = H*x1
    H = reshape([x(1:8);1], 3, 3)';
    x1 = [x(  9:9+n-1)'
          x(9+n:9+2*n-1)'];
    
    Hx1 = H*[x1
             ones(1,n)];
    x2  = Hx1(1:2,:)./Hx1(3,:);
    
    y = [x
         x2(1,:)';
         x2(2,:)'];
    
end


function [err, J] = objfun(x, m1, m2, w)
    
    n = size(m1,2);
    
    H  = reshape([x(1:8);1], 3, 3)';
    x1 = [x(  9:9+n-1)'
          x(9+n:9+2*n-1)'];
    
    Hx1 = H*[x1
             ones(1,n)];
    ih3 = 1./Hx1(3,:);
    x2  = Hx1(1:2,:) .* ih3;
    
    err1 = m1 - x1;
    err2 = m2 - x2;
    err  = [err1', err2'] .* w;

    if nargout > 1
        
        Onn = sparse(n,n);
        On8 = sparse(n,8);
        On3 = sparse(n,3);
        In1 = ones(n,1);
        Inn = spd(ones(n,1));
        x1  = x1';
        u2  = x2(1,:)';
        v2  = x2(2,:)';
        ih3 = ih3';
        u2ih3= u2.*ih3;
        v2ih3= v2.*ih3;
        de1dh = [On8
                 On8];
        de2dh = [ -ih3.*x1,      -ih3,  On3, u2ih3.*x1
                       On3,  -ih3.*x1, -ih3, v2ih3.*x1];
        de1dx = [-Inn,  Onn
                  Onn, -Inn];
        de2dx = [ spd(-H(1,1)*ih3 + H(3,1)*u2ih3),  spd(-H(1,2)*ih3 + H(3,2)*u2ih3)
                  spd(-H(2,1)*ih3 + H(3,1)*v2ih3),  spd(-H(2,2)*ih3 + H(3,2)*v2ih3) ];      
        w2 = [w; w];
        W2 = spd(w2);
        J = [de1dh,     W2*de1dx
             de2dh.*w2, W2*de2dx];
    end
       
 
end


function D = spd(v)
    if isrow(v)
        v = v';
    end
    n = length(v);
    D = spdiags(v,0,n,n);
end