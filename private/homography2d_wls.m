function H = homography2d_wls(m1, m2, normalize, w)
    
    narginchk(2,4);
    
    n = size(m1,2);
    m1 = checkpts(m1);
    m2 = checkpts(m2);
    

    if nargin < 3, normalize = true; end
    if normalize
        [m1, T1] = normalise2dpts(m1);
        [m2, T2] = normalise2dpts(m2);
    else
        T1 = eye(3);
        T2 = eye(3);
    end
    
    if nargin < 4 || isempty(w)
        w = ones(n,1);
    elseif isrow(w)
        w = w';
    end
    
    x1 = m1(1,:)';
    y1 = m1(2,:)';
    x2 = m2(1,:)';
    y2 = m2(2,:)';
    O = zeros(n,1);
    l = ones(n,1);
    
    %      H11,  H12, H13, H21,  H22, H23,      H31,     H32, H33
    A = [    O,    O,   O, -x1,  -y1,  -l,   x1.*y2,  y1.*y2,  y2
            x1,   y1,   l,   O,    O,   O,  -x1.*x2, -y1.*x2, -x2];


    M = A'*([w; w].*A);
    [~,~,V] = svd(M);
    H = reshape(V(:,9), 3, 3)';
    
    
    H = T2\H*T1;
    H = H./norm(H,'fro');    
    
end