function err = homoDist_sampson(h, m1, m2, weight)
    if nargin < 4
        weight = ones(1, size(m1,2));
    end
    if isrow(weight)
        weight = weight';
    end
    
    if length(h) == 8
        H = reshape([h;1], 3, 3)';
    else
        H = reshape(h, 3, 3)';
    end

    d   = vgg_H_sampson_distance_sqr(H, m1, m2);
    err = sqrt(d)' .* weight;
    
end