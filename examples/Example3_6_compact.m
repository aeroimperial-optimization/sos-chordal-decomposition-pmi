
% =================================================================
% Code for Example 3-6 in our paper
% Y. Zheng, G. Fantuzzi, Sum-of-squares chordal decomposition of 
%                         polynomial matrix inequalities
%
% Assume the option sos.csp in Yalmip has been modified
% =================================================================

clear;
Dim = 5:2:16;     % Dimension of the polynomial matrix
                  % the range of 10:10:80 was used in the paper

TimeSolver = zeros(length(Dim),2);
Cost       = zeros(length(Dim),2);

opts =  sdpsettings('solver','mosek');
opts.verbose = 0;
%opts.sos.newton = 0;     % no newton
%opts.sos.congruence = 0; % no symmetry

DimFull = length(Dim);
DimDec  = length(Dim);

sdpvar x y lambda

% constriants for the compact region   
g1 = 1 -x^2;
g2 = x^2 - y^2;

degree = 1;      % degree of SOS multipliers
monobasis = monolist([x,y],degree);

for index = 1:length(Dim)

    m = Dim(index);
    v = sdpvar(m,1);
    B = zeros(m);
    for i = 1:m
        for j = 1:m
            if i - j == 1 || j - i == 1
                B(i,j) = 1;
            end
        end
    end
    P = (10 + 2*x^2-x^4)*eye(m) + B*(x + x*y - x^3);

    %% Method 1: standard SOS
    % generate sos multipliers
    if index <= DimFull
        try 
            monobasis_Matrix = kron(eye(m),monobasis);
            N = length(monobasis);
            Q0 = sdpvar(m*N);
            S  = monobasis_Matrix'*Q0*monobasis_Matrix;   % SOS matrix
            s0 = v'*S*v;
            Q1 = sdpvar(m*N);
            S  = monobasis_Matrix'*Q1*monobasis_Matrix;   % SOS matrix
            s1 = v'*S*v;
            Q2 = sdpvar(m*N);
            S  = monobasis_Matrix'*Q2*monobasis_Matrix;   % SOS matrix
            s2 = v'*S*v;

            F = [sos(s0),sos(s1),sos(s2)];
            F = [F, sos(v'*(P+lambda*eye(m))*v - s0 - s1*g1 - s2*g2)];
            sosQvar = [vec(Q0);vec(Q1);vec(Q2)];
            opts.sos.csp = 0;
            sol =  solvesos(F,lambda,opts,[sosQvar;lambda]);
            TimeSolver(index,1) = sol.solvertime;
            Cost(index,1) = value(lambda);
        catch
             DimFull = index;
             warning('Out of memory');
        end
    end
    
    %% Method 2: exploiting chordal sparsity
    % generate sos multipliers
    clear s0 s1 s2 Q0 Q1 Q2
    tic    % the following construction seems very time consuming
    s0 = cell(m-1,1); Q0 = cell(m-1,1);
    s1 = cell(m-1,1); Q1 = cell(m-1,1);
    s2 = cell(m-1,1); Q2 = cell(m-1,1);

    monobasis_Matrix = kron(eye(2),monobasis);
    N = length(monobasis);
    F = [];          % SOS constraint
    temp = 0;
    sosQvar = [];     % Gram matrix Q
    for t = 1:m-1
        Q0{t} = sdpvar(2*N);
        S = monobasis_Matrix'*Q0{t}*monobasis_Matrix;   % SOS matrix
        s0{t} = v(t:t+1)'*S*v(t:t+1);

        Q1{t} = sdpvar(2*N);
        S = monobasis_Matrix'*Q1{t}*monobasis_Matrix;   % SOS matrix
        s1{t} = v(t:t+1)'*S*v(t:t+1);

        Q2{t} = sdpvar(2*N);
        S = monobasis_Matrix'*Q2{t}*monobasis_Matrix;   % SOS matrix
        s2{t} = v(t:t+1)'*S*v(t:t+1);

        F = [F,sos(s0{t}),sos(s1{t}),sos(s2{t})];
        temp = temp + s0{t} + s1{t}*g1 + s2{t}*g2;
        sosQvar = [sosQvar;vec(Q0{t});vec(Q1{t});vec(Q2{t})];
    end
    cons_time = toc;

    F = [F, sos(v'*(P+lambda*eye(m))*v - temp)];
    opts.sos.csp = 1;
    sol1 =  solvesos(F,lambda,opts,[sosQvar;lambda]);
    TimeSolver(index,2) = sol1.solvertime;
    Cost(index,2) = value(lambda);
    
    %% print information
    fprintf('m = %d,  time %5.3f  %5.3f  %5.3f\n ',m, TimeSolver(index,1),TimeSolver(index,2),cons_time);
end

