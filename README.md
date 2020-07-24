# Sum-of-squares chordal decomposition of polynomial matrix inequalities: Examples & MATLAB scripts

This repository contains MATLAB scripts that implement numerical examples for sum-of-squares (SOS) chordal decomposition results for polynomial matrices presented in the following paper:

1) Zheng, Y, and Fantuzzi, G. (2020) [ Sum-of-squares chordal decomposition of polynomial matrix inequalities ](https://arxiv.org/abs/2007.11410) arXiv preprint arXiv:2007.11410.

## Instructions

To run the scripts in this repository you need a working MATLAB installation. In addition, please install:
* The optimization toolbox [YALMIP](https://yalmip.github.io/)
* An SDP solver compatible with YALMIP. We recommend [MOSEK](https://www.mosek.com/), but any of the solvers listed [here](https://yalmip.github.io/allsolvers/) should work.

We adapted the sos.csp option in YALMIP to exploit chordal sparsity described in our paper. For this option, please
* Copy corrsparsity.m and cliquesFromSpMatD.m to the folder of /modules/sos, and replace the original corrsparsity.m (cliquesFromSpMatD.m is copied from the SparseCoLO package)

<!--- OLD CONTENT - TO BE UPDATED?
## Further notes
In the SOSexampleYalmip.m, we demonstrated three methods
* Method 1 is the standard SOS approach.
* Method 2 is based on SOS matrix decompsotion; See [2] for details.
* Method 3 is based on the correlative sparsity technique, which was orginally proposed by Waki et al 2006; See [1] for a comparsion with DSOS/SDSOS techniques.
%
Also, note that 
* Methods 2 & 3 are equivalent and are in general more conservative than Method 1 
* Methods 2 & 3 are much more scalable than Method 1 for sparse instances. See [1],[2] for more numerical examples
%
-->

## Additional references
1) Zheng, Y., Fantuzzi, G., & Papachristodoulou, A. (2018).  [ Sparse sum-of-squares (SOS) optimization: A bridge between DSOS/SDSOS and SOS optimization for sparse polynomials](https://arxiv.org/pdf/1807.05463.pdf). arXiv preprint arXiv:1807.05463.
2) Zheng, Y., Fantuzzi, G., & Papachristodoulou, A. (2018, December).  [ Decomposition and completion of sum-of-squares matrices](https://arxiv.org/pdf/1804.02711.pdf). In 2018 IEEE Conference on Decision and Control (CDC) (pp. 4026-4031). IEEE.
3) Waki, H., Kim, S., Kojima, M., & Muramatsu, M. (2006). Sums of squares and semidefinite program relaxations for polynomial optimization problems with structured sparsity. SIAM Journal on Optimization, 17(1), 218-242.

## Troubleshooting
If you have any trouble running the scripts in this repository, please email [Yang Zheng](mailto:zhengy@g.harvard.edu?Subject=SOS-csp).
