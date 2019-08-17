#                    A demo for the SOS matrix decompsotion  and Correlative Sparsity tehcnique

 Details can be found in the following papers
 [1] Zheng, Y., Fantuzzi, G., & Papachristodoulou, A. (2018).  Sparse sum-of-squares (SOS) optimization: A bridge between DSOS/SDSOS 
     and SOS optimization for sparse polynomials. arXiv preprint arXiv:1807.05463.
 [2] Zheng, Y., Fantuzzi, G., & Papachristodoulou, A. (2018, December).  Decomposition and completion of sum-of-squares matrices. 
     In 2018 IEEE Conference on Decision and Control (CDC) (pp. 4026-4031). IEEE.

% -----------------------------------------------------------
% How to use this file
% -----------------------------------------------------------

% For the correlative sparsity tehcnique, we adapted the csp option in YALMIP, 
% where we mainly changed the function sos/corrsparsity.m
% You can copy corrsparsity.m and cliquesFromSpMatD.m to the folder of
% /modules/sos, and replace the orginal corrsparsity.m


% -----------------------------------------------------------
% Before you run this file, some further instructions
% -----------------------------------------------------------

% Method 1 is the standard SOS approach:
% Method 2 is based on SOS matrix decompsotion; See [1] for details
% Method 3 is based on the correlative sparsity tehcnique
%             See [2] for a comparsion with DSOS/SDSOS techniques

% Note that 
% 1) Methods 2 & 3 are equivalent and are in general more conservative than Method 1 
% 2) Methods 2 & 3 are much more scalable than Method 1 for sparse instances
%    See [1],[2] for more numerical examples
