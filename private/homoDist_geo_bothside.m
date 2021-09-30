function [err, J] = homoDist_geo_bothside(h, m1, m2, weight)
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
    
    
    [err1, J1] = homoDist_geo_oneside(h, m1, m2, weight);
    
    if length(h) == 8
        H = reshape([h;1], 3, 3)';
    else
        H = reshape(h, 3, 3)';
    end

    
    
    adjH = [cross( H(2,:), H(3,:))
            cross( H(3,:), H(1,:))
            cross( H(1,:), H(2,:)) ]';
    Hm2 = adjH*m2;
    
    h3  = 1./Hm2(3,:);
    uv  = Hm2(1:2,:) .* h3;
    
    err2 = ( m1(1:2,:) - uv ) .* weight;
    err2 = err2';
    
    err = [err1, err2];
    
    if nargout > 1        
        dh1 = ...
            [       0,       0,       0,       0,  H(3,3), -H(3,2),       0, -H(2,3),  H(2,2)
                    0, -H(3,3),  H(3,2),       0,       0,       0,       0,  H(1,3), -H(1,2)
                    0,  H(2,3), -H(2,2),       0, -H(1,3),  H(1,2),       0,       0,       0];
        dh2 = ...
            [       0,       0,       0, -H(3,3),       0,  H(3,1),  H(2,3),       0,       0
               H(3,3),       0, -H(3,1),       0,       0,       0, -H(1,3),       0, -H(2,1)
              -H(2,3),       0,  H(2,1),  H(1,3),       0, -H(1,1),       0,       0,  H(1,1)];
        dh3 = ...
            [       0,       0,       0,  H(3,2), -H(3,1),       0, -H(2,2),  H(2,1),       0
              -H(3,2),  H(3,1),       0,       0,       0,       0,  H(1,2), -H(1,1),       0
               H(2,2), -H(2,1),       0, -H(1,2),  H(1,1),       0,       0,      0,        0];
        u = uv(1,:)';
        v = uv(2,:)';

        q = (h3.*weight)';

        J2 = ...
            [( (m2'*dh3).*u - (m2'*dh1) ).*q
             ( (m2'*dh3).*v - (m2'*dh2) ).*q];
         
        if length(h)==8
            J2 = J2(:,1:8);
        end
        
        J = [J1
             J2];
         

    end
end