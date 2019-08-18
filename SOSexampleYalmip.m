% -------------------------------------------------------------------------
%                   A demo for the SOS matrix decompsotion 
%                     and Correlative Sparsity technique
% -------------------------------------------------------------------------
% Details can be found in the following papers
% [1] Zheng, Y., Fantuzzi, G., & Papachristodoulou, A. (2018). 
%     Sparse sum-of-squares (SOS) optimization: A bridge between DSOS/SDSOS 
%     and SOS optimization for sparse polynomials. arXiv preprint arXiv:1807.05463.
% [2] Zheng, Y., Fantuzzi, G., & Papachristodoulou, A. (2018, December). 
%     Decomposition and completion of sum-of-squares matrices. 
%     In 2018 IEEE Conference on Decision and Control (CDC) (pp. 4026-4031). IEEE.

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

% This script may take around 2 minutes to run, depending on your machine,
% But you can always change the matrix dimension

clc;clear;

%Num = 10:5:40;                       % matrix dimension
Num = 20:5:40;
ResultLower = zeros(length(Num),3);   % Matrix - orignal; Matrix decomposition; Scalar-csp
ResultTime  = zeros(length(Num),3);

for Index = 1:length(Num)
    
    N = Num(Index);
    n = 2;                    % number of variables in each polynomial
    degree = 2;               % degree of each polynomial

    %% generate data
    x      = sdpvar(n,1);          % polynomial variables
    Lower  = sdpvar(1,1);          % decision variables
    mBasis = monolist(x,degree/2);

    N1  = length(mBasis);
    tmp = rand(N1);
    tmp = tmp*tmp';     % PSD matrix
    P   = mBasis'*tmp*mBasis;
    for i = 2:N
        tmp = rand(N1);
        tmp = tmp*tmp';     % PSD matrix
        P = blkdiag(P,mBasis'*tmp*mBasis);
    end

    for i = 2:N
        P(1,i) = rand(1,N1)*mBasis;
        P(i,1) = P(1,i);
    end

    %% Method 1: Matrix-original
    opts = sdpsettings('solver','mosek');
    %opts = sdpsettings('solver','sedumi');
    F = sos(P+Lower*eye(N));
    [sol,v,Q] = solvesos(F,Lower,opts);

    %% method 2: based on decomposition
    % see the paper: https://arxiv.org/pdf/1804.02711.pdf
    
    Lower1 = sdpvar(1,1);
    Pi = cell(N-1,1);
    for i = 1:N-1
        Pi{i}      = sdpvar(2);
        Pi{i}(1,2) = P(1,i+1);Pi{i}(2,1) = Pi{i}(1,2); Pi{i}(2,2) = P(i+1,i+1) + Lower1;
        Pi{i}(1,1) = polynomial(x,degree);
    end

    % constraint
    F   = [];
    tmp = 0;
    for i = 1:N-1
        tmp = tmp + Pi{i}(1,1);
        F   = [F, sos(Pi{i})];
    end
    F = [F, coefficients(tmp-P(1,1)-Lower1,mBasis) == 0];
    [sol1,v,Q] = solvesos(F,Lower1,opts);

    %% method 3: converted into scalar polynomials
     % and exploit the correlative sparsity technique
     % see the paper: https://arxiv.org/abs/1807.05463 for a comparsion with DSOS/SDSOS techniques
     
    z      = sdpvar(N,1);
    Lower3 = sdpvar(1,1);
    P3     = z'*P*z+Lower3*z'*z;
    F      = sos(P3);
    opts.sos.csp = 1;                  % only this part matters; correlative sparsity pattern
    [sol3,v3,Q3] = solvesos(F,Lower3,opts);

    ResultLower(Index,:) = [value(Lower),value(Lower1),value(Lower3)];
    ResultTime(Index,:)= [sol.solvertime,sol1.solvertime,sol3.solvertime];
end

figure;
h1 = plot(Num,ResultTime(:,1),'linewidth',2); hold on;
h2 = plot(Num,ResultTime(:,2),'linewidth',2); 
h3 = plot(Num,ResultTime(:,3),'linewidth',2); 
xlabel('Matrix dimension');ylabel('Time consumption');
legend([h1,h2,h3],'Standard SOS','SOS matrix decompostion', 'Correlative sparsity technique');
box off;
