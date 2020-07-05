% =================================================================
% Code for Example 3-5b in our paper
% Y. Zheng, G. Fantuzzi, Sum-of-squares chordal decomposition of 
%                         polynomial matrix inequalities
%
% Assume the option sos.csp in Yalmip has been modified
% =================================================================

clc;clear;
Dim = 5:5:30;     % Dimension of the polynomial matrix
                  % the range of 20:20:200 was used in the paper

TimeSolver = zeros(length(Dim),2);
Cost       = zeros(length(Dim),2);
DimFull    = length(Dim);
DimDec     = length(Dim);

sdpvar x y          % sos variables
sdpvar a b          % parameters
opts =  sdpsettings('solver','mosek');

for index = 1:length(Dim)
    m = Dim(index);
    v = sdpvar(m,1);    % vector for scalarization  

    % generating the polynomial matrix for testing
    B = zeros(m);
    for i = 1:m
        for j = 1:m
            if i - j == 1 || j - i == 1
                B(i,j) = 1;
            end
        end
    end
    P = eye(m)*x^4 + (B*(a + b) + eye(m))*x^2*y^2 + eye(m)*y^4;

    % Standard SOS without exploiting chordal sparity
    opts.sos.csp = 0;              % assume this option has been modified in Yalmip
    opts.verbose = 0;
    F = sos(v'*P*v);               % The current Yalmip only exploit csp in the scalar case 
    if index <= DimFull
        try
            sol = solvesos(F,a+b,opts);
            TimeSolver(index,1) = sol.solvertime;
            Cost(index,1)       = value(a+b);
        catch
             DimFull = index;
             warning('Out of memory');
        end
    end

     % SOS via exploiting chordal sparity
     if index <= DimDec
        opts.sos.csp = 1;              % only this option matters
        F = sos(v'*P*v);               % The current Yalmip only exploit csp in the scalar case     
        try
            sol1 = solvesos(F,a+b,opts);
            TimeSolver(index,2) = sol1.solvertime;
            Cost(index,2) = value(a+b);
        catch
            DimDec = index;
            warning('Out of memory');
        end
    end
    
    fprintf('m = %d,  time %5.3f  %5.3f \n ',m, TimeSolver(index,1),TimeSolver(index,2));
end

