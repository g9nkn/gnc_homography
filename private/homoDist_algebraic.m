function [err, J] = homoDist_algebraic(h, m1, m2, weight)
    if nargin < 4
        weight = ones(1, size(m1,2));
    end
    if isrow(weight)
        weight = weight';
    end
    
    if length(h) == 8
        h = [h; 1];
    end
%     h = h/norm(h);

    
    x1 = m1(1,:)';
    y1 = m1(2,:)';
    x2 = m2(1,:)';
    y2 = m2(2,:)';
    
    n = size(m1,2);
    O = zeros(n, 1);
    l = ones(n, 1);
     
    %       H11,  H12, H13, H21,  H22, H23,     H31,     H32, H33,
    M1 = [     O,    O,   O, -x1,  -y1,  -l,  x1.*y2,  y1.*y2,  y2];
    M2 = [    x1,   y1,   l,   O,    O,   O, -x1.*x2, -y1.*x2, -x2];
    err = [M1*h, M2*h] .* weight;
    
    if nargout > 1
        %J = autoJac(@homoDist_sampson, h, m1, m2, weight); % doesn't work with reshape
    end
    
%     h = [h
%          1];
%     n = size(m1,2);
%     
%     x1 = m1(1,:)';
%     y1 = m1(2,:)';
%     x2 = m2(1,:)';
%     y2 = m2(2,:)';
%     
%     O = zeros(n, 1);
%     l = ones(n, 1);
%      
%     %       H11,  H12, H13, H21,  H22, H23,     H31,     H32, H33,
%     M1 = [     O,    O,   O, -x1,  -y1,  -l,  x1.*y2,  y1.*y2,  y2];
%     M2 = [    x1,   y1,   l,   O,    O,   O, -x1.*x2, -y1.*x2, -x2];
%                   
%     N11 = [     O,    O,   O,  -l,    O,   O,      y2,       O,   O];
%     N12 = [     O,    O,   O,   O,   -l,   O,       O,      y2,   O];
%     N13 = [     O,    O,   O,   O,    O,   O,      x1,      y1,   l];
%     N21 = [     l,    O,   O,   O,    O,   O,     -x2,       O,   O];
%     N22 = [     O,    l,   O,   O,    O,   O,       O,     -x2,   O];
%     N23 = [     O,    O,   O,   O,    O,   O,     -x1,     -y1,  -l];
%     
%     alg1 = M1*h;
%     alg2 = M2*h;
%     den1 = sqrt( (N11*h).^2 + (N12*h).^2 + (N13*h).^2 );
%     den2 = sqrt( (N21*h).^2 + (N22*h).^2 + (N23*h).^2 );
%     
%     
%     err = [alg1./den1, alg2./den2] .* weight;
%     
%     if nargout > 1
%     end
end