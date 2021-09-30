function [err, J] = homoDist_geo_oneside(h, m1, m2, weight)
    
    if size(m1,1) == 2
        m1 = [m1
              ones(1,size(m1,2))];
    end
    if size(m2,1) == 2
        m2 = [m2
              ones(1,size(m2,2))];
    end    
    
    if nargin < 4
        weight = ones(1, size(m1,2));
    end
    if iscolumn(weight)
        weight = weight';
    end
    
    
    if length(h) == 8
        H = reshape([h;1], 3, 3)';
    else
        H = reshape(h, 3, 3)';
    end

    
    Hm1 = H*m1;
    h3  = 1./Hm1(3,:);
    uv  = Hm1(1:2,:) .* h3;
    
    err = ( m2(1:2,:) - uv ) .* weight;
    err = err';
    
    if nargout > 1
        u  = uv(1,:)';
        v  = uv(2,:)';
        s  = (h3 .* weight)';
        
        n = size(m1,2);
        O = zeros(n,3);
            
        Q = -m1' .* s;
        J = [Q, O, -Q.*u
             O, Q, -Q.*v];
        
        if length(h) == 8
            J = J(:,1:8);
        end
        
    end
end