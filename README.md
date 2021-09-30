# Homography estimation using adaptive GNC

&copy; 2021 NEC Corporation

This repository is an official MATLAB implementation of the paper "[Fast and Robust Homography Estimation by Adaptive Graduated Non-Convexity", ICIP2019](https://jpn.nec.com/rd/people/docs/icip2019_nakano_shibata.pdf) [\[1\]](#reference).

## License

This software is released under the NEC Corporation License.
See [LICENSE](https://github.com/g9nkn/gnc_homography/LICENSE) before using the code. If you use this code, please cite the paper.

```bibtex
@inproceedings{nakano2019fast,
  title={Fast and robust homography estimation by adaptive graduated non-convexity},
  author={Nakano, Gaku and Shibata, Takashi},
  booktitle={2019 IEEE International Conference on Image Processing (ICIP)},
  pages={2556--2560},
  year={2019},
  organization={IEEE}
}
```

For commercial use, please contact Gaku Nakano \<g-nakano@nec.com\>.

## Usage

### Install

Copy `gnc_homography` folder and set a path to folders by `addpath('gnc_homography','gnc_homography/3rdparty')`.

### Function API

```
[H, inliers] = gnc_homog_nakano_icip2019(m1, m2, 'param1', value1, 'param2', value2, ...)

INPUTS
  m1, m2     : (2xN or 3xN matrix)
               2D point correspondences in inhomogeneous coordinates
               [x1, x2, ... xN
                y1, y2, ... yN]
                  or
               homogeneous coordinates (will be normalized by [xi/zi; yi/zi])
               [x1, x2, ... xN
                y1, y2, ... yN
                z1, z2, ... zN]

OPTIONS
  lambda_max : (scalar, default: 1e3)
               The maxilambdam pixel error for determining inliers.

  lambda_min : (scalar, default: 2)
               The minilambdam pixel error for determining inliers.

  rho        : (scalar, default: 0.95)
               Decreasing factor for updating lambda = rho*lambda in each iteration.
               If lambda_new - lambda < 0.5 for any rho, lambda_new = lambda - 0.5.

  method     : (char, default: 'symtrans')
               Homography estimation method which minimizes
               'symtrans'      : symmetric transfer error
               'singletrans'   : error in the first image
               'reproj'        : reprojection error
               'sampson'       : Sampson error (1st approximation of reprojection error)
               'algebraic'     : algebraic error

  weightfunc : (char, default: 'truncL2') 
               Weighting function of the IRLS scheme
               'truncL2'       : truncated L2
               'Tukey'         : Tukey's biweight
               'GemanMcClure'  : Geman-McClure
               'L2log'         : L2-log
               'Entropy'       :
               'MFT'           : Mean-Field function

  w          : (Nx1 or 1xN vector, defalut: ones(n,1)) 
               Initial weight for the point correspondences (0<=wi<=1).
               Given this parameter, the initial homography H is computed by weighted least squares and
               lambda_max is overwritten based on residuals of w > 0.5 points. Ohterwise, H is given by eye(3).

  choosebest : (logical, default: true)
               Return the best model if true. Otherwise, return the
               model at the final GNC iteration.

  verbose    : (logical, default: true)
               Display outputs at each GNC itereation. 

OUTPUTS
  H          : (3x3 matrix)
               Homography matrix such that m2 = H*m1.

  inliers    : (1xN logical vector)
               Indeces of inlier points.
```

Run `demo/demo_gnc_homog_icip2019.m` to understand how it works.  
Please refer to [\[2\]](#reference) and [\[3\]](#reference) for the details of GNC, IRLS (M-estimator), and weight functions.

## Reference

1. Gaku Nakano and Shibata Takashi, "Fast and Robust Homography Estimation by Adaptive Graduated Non-Convexity," ICIP 2019. [(pdf)](https://jpn.nec.com/rd/people/docs/icip2019_nakano_shibata.pdf)

2. Michael J. Black and Anand Rangarajan, "On the unification of line processes, outlier rejection, and robust statistics with applications in early vision," IJCV, 19(1):57-91, 1996, Springer. [(pdf)](https://www.cise.ufl.edu/~anand/pdf/ijcv.pdf)

3. Zhengyou Zhang, "Parameter estimation techniques: A tutorial with application to conic fitting," IVC, 15(1):59-76, 1997, Elsevier. [(pdf)](https://www.microsoft.com/en-us/research/wp-content/uploads/2016/11/RR-2676.pdf)

## Contributors

- Gaku Nakano, Central Research Laboratories, NEC Corporation.  
<g-nakano@nec.com>  
(ENG) <https://www.nec.com/en/global/rd/people/gaku_nakano.html>  
(JPN) <https://jpn.nec.com/rd/people/gaku_nakano.html>

- Takashi Shibata, Central Research Laboratories, NEC Corporation.  
<t.shibata@ieee.org>
