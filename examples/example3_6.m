% =================================================================
% Code for Example 3-6 in the paper:
% Y. Zheng, G. Fantuzzi, Sum-of-squares chordal decomposition of polynomial 
% matrix inequalities
%
% Assume the option sos.csp in Yalmip has been modified
% =================================================================

% Clean up
clear
yalmip clear

% Parameters
Deg = 1;             % degree of SOS multipliers
                     % the range of 1:4 was used in the paper
Dim = 5:5:15;        % Dimension of the polynomial matrix
                     % the range of 10:10:80 was used in the paper
  
% Initialize containers
TimeSolver = cell(2,1);
Cost       = cell(2,1);
TimeSolver{1} = zeros(length(Dim),2);
Cost{1}       = zeros(length(Dim),2);
TimeSolver{2} = zeros(length(Dim),2);
Cost{2}       = zeros(length(Dim),2);

% YALMIP options and variables
opts =  sdpsettings();
opts.verbose = 0;
sdpvar x y lambda lambda1

% Loop over degrees
for ind = 1:length(Deg)
    TimeSolver{ind} = zeros(length(Dim),2);
    Cost{ind} = zeros(length(Dim),2);
    DimFull = length(Dim);
    DimDec = length(Dim);

    % constriants for the compact region   
    g1 = 1 -x^2;
    g2 = x^2 - y^2;
    degree = Deg(ind);
    monobasis = monolist([x,y],degree);

    for index = 1:length(Dim)

        m = Dim(index);
        v = sdpvar(m,1);     % scalarization
        B = zeros(m);
        B(:,1) = 1; B(1,:) = 1; B(1,1) = 0;
        P = (10 + y^3-x^4)*eye(m) + B*(x + x*y - x^3);

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
                F = [F, sos(v'*(P-lambda*eye(m))*v - s0 - s1*g1 - s2*g2)];
                sosQvar = [vec(Q0);vec(Q1);vec(Q2)];
                opts.sos.csp = 0;
                sol =  solvesos(F,-lambda,opts,[sosQvar;lambda]);
                TimeSolver{ind}(index,1) = sol.solvertime;
                Cost{ind}(index,1) = value(lambda);
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
        F = [];           % SOS constraint
        temp = 0;
        sosQvar = [];     % Gram matrix Q
        for t = 1:m-1
            Q0{t} = sdpvar(2*N);
            S = monobasis_Matrix'*Q0{t}*monobasis_Matrix;   % SOS matrix
            s0{t} = v([1,t+1],1)'*S*v([1,t+1],1);

            Q1{t} = sdpvar(2*N);
            S = monobasis_Matrix'*Q1{t}*monobasis_Matrix;   % SOS matrix
            s1{t} = v([1,t+1],1)'*S*v([1,t+1],1);

            Q2{t} = sdpvar(2*N);
            S = monobasis_Matrix'*Q2{t}*monobasis_Matrix;   % SOS matrix
            s2{t} = v([1,t+1],1)'*S*v([1,t+1],1);

            F = [F,sos(s0{t}),sos(s1{t}),sos(s2{t})];
            temp = temp + s0{t} + s1{t}*g1 + s2{t}*g2;
            sosQvar = [sosQvar;vec(Q0{t});vec(Q1{t});vec(Q2{t})];
        end
        cons_time = toc;

        F = [F, sos(v'*(P-lambda1*eye(m))*v - temp)];
        opts.sos.csp = 1;
        sol1 =  solvesos(F,-lambda1,opts,[sosQvar;lambda1]);
        TimeSolver{ind}(index,2) = sol1.solvertime;
        Cost{ind}(index,2) = value(lambda1);

        %% print information
        fprintf('m = %d,  time %5.3f  %5.3f  %5.3f\n ',m, TimeSolver{ind}(index,1),TimeSolver{ind}(index,2),cons_time);
    end
end

