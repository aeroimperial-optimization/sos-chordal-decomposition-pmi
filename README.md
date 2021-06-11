# Sum-of-squares chordal decomposition of polynomial matrix inequalities: Examples

This repository contains MATLAB scripts that implement the numerical examples in the following paper:

1) Zheng, Y. and Fantuzzi, G. (2020). Sum-of-squares chordal decomposition of polynomial matrix inequalities, [arXiv:2007.11410 [math.OC]](https://arxiv.org/abs/2007.11410).

## Instructions

To run the scripts in this repository you need a working MATLAB installation. In addition, please install:
* The optimization toolbox [YALMIP](https://yalmip.github.io/)
* An SDP solver compatible with YALMIP. We recommend [MOSEK](https://www.mosek.com/), but any of the solvers listed [here](https://yalmip.github.io/allsolvers/) should work.

We adapted the sos.csp option in YALMIP to exploit chordal sparsity described in our paper. For this option, please
* Copy corrsparsity.m and cliquesFromSpMatD.m to the folder of /modules/sos, and replace the original corrsparsity.m (cliquesFromSpMatD.m is copied from the SparseCoLO package)

## Additional references
1) Zheng, Y., Fantuzzi, G. and Papachristodoulou, A. (2019). Sparse sum-of-squares (SOS) optimization: A bridge between DSOS/SDSOS and SOS optimization for sparse polynomials. In: *Proceedings of the 2019 American Control Conference*, pp. 5513-5518.
[[Publisher Link](https://doi.org/10.23919/ACC.2019.8814998)] [[Open Access Link](https://arxiv.org/pdf/1807.05463.pdf)]

2) Zheng, Y., Fantuzzi, G. and Papachristodoulou, A. (2018). Decomposition and completion of sum-of-squares matrices. In: *Proceedings of the 57th IEEE Conference on Decision and Control*, pp. 4026-4031.
[[Publisher Link](https://doi.org/10.1109/CDC.2018.8619144)] [[Open Access Link](https://arxiv.org/pdf/1804.02711.pdf)]

3) Waki, H., Kim, S., Kojima, M. and Muramatsu, M. (2006). Sums of squares and semidefinite program relaxations for polynomial optimization problems with structured sparsity. *SIAM Journal on Optimization*, **17**(1), 218-242.
[[Publisher Link](https://doi.org/10.1137/050623802)]

## Troubleshooting
If you have any trouble running the scripts in this repository, please email [Yang Zheng](mailto:zhengy@g.harvard.edu?Subject=SOS-csp).
