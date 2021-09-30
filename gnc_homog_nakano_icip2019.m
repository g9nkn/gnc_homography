%% A MATLAB code of the homography estimation using adaptive GNC [1].
% Copyright (c) 2021 NEC Corporation
% This software is released under the NEC Corporation License. See the License file.
% For commercial use, please contact Gaku Nakano <g-nakano@nec.com>.
%
% USAGE
%   [H, inliers, stats] = gnc_homog_nakano_icip2019(m1, m2, 'param1', value1, 'param2', value2, ...)
%
% INPUTS
%   m1, m2     : (2xN or 3xN matrix)
%                2D point correspondences in inhomogeneous coordinates
%                [x1, x2, ... xN
%                 y1, y2, ... yN]
%                   or
%                homogeneous coordinates (will be normalized by [xi/zi; yi/zi])
%                [x1, x2, ... xN
%                 y1, y2, ... yN
%                 z1, z2, ... zN]
%
% OPTION PARAMETERS
%   lambda_max : (scalar, default: 1e3)
%                The maxilambdam pixel error for determining inliers.
%
%   lambda_min : (scalar, default: 2)
%                The minilambdam pixel error for determining inliers.
%
%   rho        : (scalar, default: 0.95)
%                Decreasing factor for updating lambda = rho*lambda in each iteration.
%                If lambda_new - lambda < 0.5 for any rho, lambda_new = lambda - 0.5.
%
%   method     : (char, default: 'symtrans')
%                Homography estimation method which minimizes
%                'symtrans'      : symmetric transfer error
%                'singletrans'   : error in the first image
%                'reproj'        : reprojection error
%                'sampson'       : Sampson error (1st approximation of reprojection error)
%                'algebraic'     : algebraic error
%
%   weightfunc : (char, default: 'truncL2') 
%                Weighting function of the IRLS scheme
%                'truncL2'       : truncated L2
%                'Tukey'         : Tukey's biweight
%                'GemanMcClure'  : Geman-McClure
%                'L2log'         : L2-log
%                'Entropy'       :
%                'MFT'           : Mean-Field function
%
%   w          : (Nx1 or 1xN vector, defalut: ones(n,1)) 
%                Initial weight for the point correspondences (0<=wi<=1).
%                Given this parameter, the initial homography H is computed by weighted least squares and
%                lambda_max is overwritten based on residuals of w > 0.5 points. Ohterwise, H is given by eye(3).
%
%   choosebest : (logical, default: true)
%                Return the best model if true. Otherwise, return the
%                model at the final GNC iteration.
%
%   verbose    : (logical, default: true)
%                Display outputs at each GNC itereation. 
%
% OUTPUTS
%   H          : (3x3 matrix)
%                Homography matrix such that m2 = H*m1.
%
%   inliers    : (1xN logical vector)
%                Indeces of inlier points.
%
%   stats      : (struct, optional)
%                
% REFERENCE:
%   [1] Gaku Nakano and Shibata Takashi, 
%       "Fast and Robust Homography Estimation by Adaptive Graduated Non-Convexity,"
%       pp. 2556-2560, ICIP 2019.  
function [H, inliers, stats, beststats] = gnc_homog_nakano_icip2019(m1, m2, varargin)
    

    %% check argments
    assert(nargin >= 2, 'Not enough argments.');

    m1 = checkpts(m1);
    m2 = checkpts(m2);
    assert( size(m1,2)==size(m2,2), 'Point size does not match.' );
    
	numPts = size(m1,2);
    opts   = getOptions(numPts, varargin{:});
    
    
    
    
    %% Given initial weight but initial H, calc initial guess of H
    w = opts.w;
    h = opts.h8x1;
    if opts.initw && ~opts.initH
        h = warpper_algsolver(m1, m2, true, w);
    end
    
    
    
    %% set IRLS parameters
    [opts, x] = setOptimizationOptions(h, m1, m2, opts);    
    lambda    = opts.lambda_max;
    maxitr    = 2 * round( (log(opts.lambda_min) - log(opts.lambda_max))/log(opts.rho) );
    itr       = 0;
    stop      = 0;
    totalErr  = inf;
    p         = w>0.5;
    if opts.initw || opts.initH
        err    = opts.errfunc(x, m1, m2);
        err2   = sum(err.^2, 2); 
        lambda = round( max(min(opts.lambda_max, 2*std2(err(p,:))+mean2(err(p,:))), opts.lambda_min) );
    else
        err2 = inf(numPts, 1);
    end

    
    
    
    %% start GNC
    if opts.stats
        stats    = createStats(numPts, maxitr+1);
        stats(1) = recordStats(w, lambda, totalErr, inf(numPts,1), inf, inf, [h;1] );
    end
    
    if opts.verbose
        fprintf('%4s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'itr', 'w.*f', 'norm(dh)', 'lambda', 'norm(dw)', 'InRatio');
        fprintf('%4d\t%10.4g\t%10.4g\t%10.4g\t%10.4g\t%10.4g (%d)\n', itr, sum(w.*sqrt(err2))/sum(w), 0, lambda, 0, nnz(p)/numPts*100, nnz(p));
    end
    
    
    while ~stop
        itr    = itr + 1;
        h_old  = h;
        w_old  = w;

        
        % IRLS scheme
        err       = opts.errfunc(x, m1, m2);               % 1) calc residuals on current model param
        err2      = sqrt(sum(err.^2, 2)) / size(err,2);                          
        w         = opts.weightfunc(err2, lambda);         % 2) determine weight
        x         = opts.optimfunc(x, sqrt(w));            % 3) optimize model param with new weight
        h         = x(1:8);
        
        
        dh       = norm(h - h_old)/norm(h_old);
        dw       = norm(w - w_old)/norm(w_old);
        p        = w > 0.5;
        if opts.stats
            stats(itr+1) = recordStats(w, lambda, totalErr, err2, dh, dw, [h;1] );
        end
        
        if opts.verbose
            fprintf('%4d\t%10.4g\t%10.4g\t%10.4g\t%10.4g\t%10.4g (%d)\n', itr, sum(w.*sqrt(err2))/sum(w), dh, lambda, dw, nnz(p)/numPts*100, nnz(p));
        end

        
        % check stopping criteria
        if lambda <= opts.lambda_min %&& dh < opts.tol
            stop = 1;
        elseif itr == maxitr
            stop = 10;
        end
        
        
        
        % update lambda
        if opts.fastconv
            err_in     = err2(p,:);
            lambda_new = round( min(opts.rho*lambda, 2*std2(err_in)+mean2(err_in) ) );
        else
            lambda_new = opts.rho*lambda;
        end
        
        if opts.fastconv && lambda - lambda_new < 0.5
            lambda_new = lambda - 0.5;
        end
        lambda = max(lambda_new, opts.lambda_min);
    end
    
    
    
    if opts.stats && itr < maxitr+1
        stats(itr+2:end) = [];
    end
    
    
    if opts.choosebest
        [H, inliers, beststats, best_iter] = selectbest(stats);
    else
        inliers = p';
        H = reshape([h;1], 3, 3)';
        best_iter = itr;
    end
    H = H./norm(H,'fro');
    
    if opts.verbose
        fprintf('Output homography was obtained at %d-th iteration\n', best_iter);
    end
    
end


function opts = getOptions(npts, varargin)
    funcList   = {'truncL2', 'GemanMcClure', 'L2log', 'Entropy',  'MFT', 'Tukey'};
    methodList = {'symtrans', 'singletrans', 'reproj', 'sampson', 'algebraic'};
    
    maxlambda  = 1e3;
    minlambda  = 2;
    rho        = 0.95;
    weightfunc = funcList{1};
    method     = methodList{1};
    weight     = ones(npts,1);
    H          = eye(3);
    verbose    = true;
    fast       = false;
    tol        = 1e-5;
    
    p = inputParser;
    addParameter(p, 'lambda_max', maxlambda , @(x)assert( isnumeric(x) && isscalar(x) && x > 0         , 'lambda_max must be a scalar greater than 0. (default:1e3)' ) );
    addParameter(p, 'lambda_min', minlambda , @(x)assert( isnumeric(x) && isscalar(x) && x > 0         , 'lambda_min must be a scalar greater than 0. (default:2)' ) );
    addParameter(p, 'rho'       , rho       , @(x)assert( isnumeric(x) && isscalar(x) && x > 0 && x < 1, 'rho must be a scalar greater than 0 and less than 1. (default: 0.95)' ) );
    addParameter(p, 'weightfunc', weightfunc, @(x)assert( ischar(x) && any(strcmp(funcList, x))        , ['weightfunc must be one of ' cell2char(funcList) '. (default: ''truncL2'')']) );
    addParameter(p, 'method'    , method    , @(x)assert( ischar(x) && any(strcmp(methodList, x))      , ['method must be one of ' cell2char(methodList) '. (default: ''symtrans'')']) );
    addParameter(p, 'H'         , H         , @(x)assert( isnumeric(x) && all([3,3]==size(x))          , 'H must be 3x3 matrix (default: eye(3))') );
    addParameter(p, 'w'         , weight    , @(x)assert( isvector(x) && numel(x)==npts && all(x>=0 & x<=1), 'w must be a Nx1 vector having 0<=w(i)<=1. (default: ones(N,1))') );
    addParameter(p, 'verbose'   , verbose   , @(x)assert( isscalar(x) && islogical(x)                  , 'verbose must be true or false. (default: true)') );
    addParameter(p, 'fastconv'  , fast      , @(x)assert( isscalar(x)                                  , 'fastconv must be true or false. (default: false)') );
    addParameter(p, 'tol'       , tol       , @(x)assert( isscalar(x) && isnumeric(x) && x > 0         , 'tol must be a scalar greater than 0. (default: 1e-5)') );
    addParameter(p, 'choosebest', true      , @(x)assert( isscalar(x) && islogical(x)                  , 'choosebest must be true or false. (default: true)') );

    parse(p, varargin{:});
    opts = p.Results;
    
    opts.weightfunc = str2func(['Lossfun_', opts.weightfunc]);
    opts.stats = true;
    if isrow(opts.w)
        opts.w = opts.w';
    end
    
    if islogical(opts.w)
        opts.w = double(opts.w);
    end
    
    if any(strcmp(varargin, 'w'))
        opts.initw = true;
    else
        opts.initw = false;
    end
    
    if any(strcmp(varargin, 'H'))
        opts.initH = true;
    else
        opts.initH = false;
    end
    opts.h8x1 = [opts.H(1,:), opts.H(2,:), opts.H(3,1:2)]' / opts.H(3,3);
    
end

function str = cell2char(c)
    str = '{';
    for i=1:length(c)
        str = [str, '''', c{i}, ''','];
    end
    str(end) = '}';
end

function [opts, x] = setOptimizationOptions(h, m1, m2, opts)
    
    opts_optim = optimoptions('lsqnonlin', 'algo', 'levenberg-marquardt', 'disp', 'none', 'SpecifyObjectiveGradient', true, 'UseParallel', false);
    switch opts.method
        case 'algebraic'
            opts.errfunc   = @homoDist_geo_oneside;
            opts.optimfunc = @(x0, w) warpper_algsolver(m1, m2, true, w);

        case 'singletrans'
            opts.errfunc   = @homoDist_geo_oneside;
            opts.optimfunc = @(x0, w) lsqnonlin(@(x)homoDist_geo_oneside(x,m1,m2,w), x0, [], [], opts_optim);

        case 'symtrans'
            opts.errfunc   = @homoDist_geo_bothside;
            opts.optimfunc = @(x0, w) lsqnonlin(@(x)homoDist_geo_bothside(x,m1,m2,w), x0, [], [], opts_optim);
            
        case 'sampson'
            opts_optim.SpecifyObjectiveGradient = false;
            opts.errfunc   = @homoDist_geo_oneside;
            opts.optimfunc = @(x0, w) lsqnonlin(@(x)homoDist_sampson(x,m1,m2,w), x0, [], [], opts_optim);
            
        case 'reproj'
            opts.errfunc   = @homoDist_reproj;
            opts.optimfunc = @(x0, w) H_reproj_lsqnonlin(x0, m1, m2, w, opts_optim);
            
        otherwise
            error('wrong method');
    end
    
    
    if strcmp(opts.method, 'reproj')
        [m1, m2] = H_triang(h, m1, m2);
        x        = [h; m1(1,:)'; m1(2,:)'; m2(1,:)'; m2(2,:)'];
    else
        x        = h;
    end
    
end



function h = warpper_algsolver(m1, m2, normalise, w)
    H = homography2d_wls(m1, m2, normalise, w);
    h = [H(1,:), H(2,:), H(3,1:2)]'/H(3,3);
end