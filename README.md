#                    A demo for the SOS matrix decompsotion  and Correlative Sparsity technique

 Details can be found in the following papers
 
 * [1] Zheng, Y., Fantuzzi, G., & Papachristodoulou, A. (2018).  [ Sparse sum-of-squares (SOS) optimization: A bridge between DSOS/SDSOS 
     and SOS optimization for sparse polynomials](https://arxiv.org/pdf/1807.05463.pdf). arXiv preprint arXiv:1807.05463.
     
* [2] Zheng, Y., Fantuzzi, G., & Papachristodoulou, A. (2018, December).  [ Decomposition and completion of sum-of-squares matrices](https://arxiv.org/pdf/1804.02711.pdf). In 2018 IEEE Conference on Decision and Control (CDC) (pp. 4026-4031). IEEE.

## How to use this file

For the correlative sparsity technique, we adapted the csp option in YALMIP, where we mainly changed the function sos/corrsparsity.m.
You can copy corrsparsity.m and cliquesFromSpMatD.m to the folder of /modules/sos, and replace the original corrsparsity.m.

* Note that cliquesFromSpMatD.m is copied from the SparseCoLO package.

## Some further notes
In the SOSexampleYalmip.m, we demonstrated three methods
* Method 1 is the standard SOS approach.
* Method 2 is based on SOS matrix decompsotion; See [2] for details.
* Method 3 is based on the correlative sparsity technique, which was orginally proposed by Waki et al 2006; See [1] for a comparsion with DSOS/SDSOS techniques.

Also, note that 
* 1) Methods 2 & 3 are equivalent and are in general more conservative than Method 1 
* 2) Methods 2 & 3 are much more scalable than Method 1 for sparse instances. See [1],[2] for more numerical examples


## Contact us<a name="Contacts"></a>
If you couldn't see computational time improvements or have any trouble running the file, please email [Yang Zheng](mailto:zhengy@g.harvard.edu?Subject=SOS-csp).
