function [x1, x2] = H_triang(h, m1, m2)
    
    n = size(m1,2);
    x1 = zeros(2,n);
    x2 = zeros(2,n);
    
    if length(h) == 8
        H = reshape([h;1], 3, 3)';
    else
        H = reshape(h, 3, 3)';
    end
    H = H/norm(H,'fro');
    
    options = optimoptions('fmincon', 'disp', 'none', 'algo', 'sqp', 'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true);
    for i = 1:n
        m  = [m1(1:2,i); m2(1:2,i)];
        x0 = [m1(1:2,i); m2(1:2,i)];
        x  = fmincon(@(x) objfunc(x, m), x0, [], [], [], [], [], [], @(x)nonlcon(x,H), options );
        x1(:,i) = x(1:2);
        x2(:,i) = x(3:4);
    end

     
end


function [err, J] = objfunc(x, m)
    err = norm(x - m)^2;
    
    if nargout > 1
    J   = 2*(x - m);
    end
end

function [c, ceq, grad_c, grad_ceq] = nonlcon(x, H)
    c      = [];
    grad_c = [];
    
    eqs = skew([x(3:4);1])*H*[x(1:2);1];
    ceq = eqs(1:2);

    if nargout > 2
        grad_ceq = [   H(3,1)*x(4) - H(2,1),   H(3,2)*x(4) - H(2,2),                              0,   H(3,3) + H(3,1)*x(1) + H(3,2)*x(2)
                       H(1,1) - H(3,1)*x(3),   H(1,2) - H(3,2)*x(3), - H(3,3) - H(3,1)*x(1) - H(3,2)*x(2),                             0]';
    end
end

function S = skew(a)
    S = [  0   -a(3)  a(2)
          a(3)   0   -a(1)
         -a(2)  a(1)   0];
end