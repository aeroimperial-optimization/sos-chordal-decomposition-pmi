% =================================================================
% Code for Example 5.1 in the paper:
% Y. Zheng, G. Fantuzzi, Sum-of-squares chordal decomposition of polynomial 
% matrix inequalities
%
% NOTE: Assume the option sos.csp in Yalmip has been modified or
%          install the following fork of YALMIP:
% https://github.com/aeroimperial-optimization/aeroimperial-yalmip
% =================================================================
clc
clear

% YALMIP variables and options
sdpvar x y z          % sos variables
sdpvar k1 k2          % parameters
opts =  sdpsettings();
opts.verbose = 0;

% Parameters
Dim = 5:5:40;     % Dimension of the polynomial matrix
                 % the range of 5:5:40 was used in the paper                  
nu = 0:1:3;          % degree of SOS multipler  
                 % the range of 1:4 was used in the paper            

% Initalize containers
TimeSolver = cell(length(nu),1);  % time consumption by Mosek
Cost       = cell(length(nu),1);  % Cost value

% Loop over nu
for dind = 1:length(nu)
    TimeSolver{dind} = zeros(length(Dim),2);
    Cost{dind} = zeros(length(Dim),2);
    DimFull = length(Dim);
    DimDec = length(Dim);

    for index = 1:length(Dim)
        d = Dim(index);

        % generating the polynomial matrix for testing
        k = [k1;k2];
        m = 3*d;
        a = [k2*x^4+y^4; k2*y^4+z^4; k2*z^4+x^4];
        b = [x^2*y^2; y^2*z^2; z^2*x^2];
        clear P;
        for i = 1:m-1
            if mod(i,2)
                c = k1;
            else
                c = k2;
            end
            ind = rem([i-1,i],3)+1;
            P(i,i) = a(ind(1));
            P(i+1,i+1) = a(ind(2));
            P(i,i+1) = b(ind(1))*c;
            P(i+1,i) = b(ind(1))*c;
        end
        clear c
        v = sdpvar(m,1);    % vector for scalarization  

        % Standard SOS without exploiting chordal sparity
        opts.sos.csp = 0;
        F = sos(v'*P*v*(x^2+y^2+z^2)^(nu(dind)));
        if index <=  DimFull
            try
                sol = solvesos(F,k2 - 10*k1,opts);
                TimeSolver{dind}(index,1) = sol.solvertime;
                Cost{dind}(index,1)       = value(k2 - 10*k1);
            catch
                 DimFull = index;
                 warning('Out of memory');
            end
        end

         % SOS via exploiting chordal sparity
         if index <= DimDec
            opts.sos.csp = 1;
            F = sos(v'*P*v*(x^2+y^2+z^2)^(nu(dind)));
            try
                sol1 = solvesos(F,k2 - 10*k1,opts);
                TimeSolver{dind}(index,2) = sol1.solvertime;
                Cost{dind}(index,2) = value(k2 - 10*k1);
            catch
                DimDec = index;
                warning('Out of memory');
            end
        end

        fprintf('m = %d,  time %5.3f  %5.3f \n',m, TimeSolver{dind}(index,1),TimeSolver{dind}(index,2));
    end
end
