% =================================================================
% Code for Example 3-5b in our paper
% Y. Zheng, G. Fantuzzi, Sum-of-squares chordal decomposition of 
%                         polynomial matrix inequalities
%
% Assume the option sos.csp in Yalmip has been modified
% =================================================================
warning off;
clc;clear;
sdpvar x y z          % sos variables
sdpvar k1 k2          % parameters
opts =  sdpsettings('solver','mosek');

Dim = 5:5:40;     % Dimension of the polynomial matrix
                  % the range of 20:20:200 was used in the paper                  
nu = 0:3;         % degree of SOS multipler  

TimeSolver = cell(length(nu),1);
Cost       = cell(length(nu),1);

for dind = 1:length(nu)
    % each degree of SOS multipler 
    
    TimeSolver{dind} = zeros(length(Dim),2);
    Cost{dind}       = zeros(length(Dim),2);
    DimFull    = length(Dim);
    DimDec     = length(Dim);

    
    for index = 1:length(Dim)
        d = Dim(index);

        % generating the polynomial matrix for testing

        k = [k1;k2];
        m = 3*d;
        a = [k2*x^4+y^4; k2*y^4+z^4; k2*z^4+x^4];
        b = [x^2*y^2; y^2*z^2; z^2*x^2];
        for i = 1:m-1
            if mod(i,2)
                c = k1;
            else
                c = k2;
            end
            ind = rem([i-1,i],3)+1;
            P(i,i) = a(ind(1));
            P(i+1,i+1) = a(ind(2));
            P(i,i+1) = b(ind(1))*c; % +temp; % dirty trick to force dense couplings in x,y,z
            P(i+1,i) = b(ind(1))*c; % +temp; % dirty trick to force dense couplings in x,y,z
        end
        clear c
        v = sdpvar(m,1);    % vector for scalarization  

        % Standard SOS without exploiting chordal sparity
        opts.sos.csp = 0;              % assume this option has been modified in Yalmip
        opts.verbose = 0;
        F = sos(v'*P*v*(x^2+y^2+z^2)^(nu(dind)));               % The current Yalmip only exploit csp in the scalar case 
        if index <= DimFull
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
            opts.sos.csp = 1;              % only this option matters
            F = sos(v'*P*v*(x^2+y^2+z^2)^(nu(dind)));               % The current Yalmip only exploit csp in the scalar case     
            try
                sol1 = solvesos(F,k2 - 10*k1,opts);
                TimeSolver{dind}(index,2) = sol1.solvertime;
                Cost{dind}(index,2) = value(k2 - 10*k1);
            catch
                DimDec = index;
                warning('Out of memory');
            end
        end

        fprintf('m = %d,  time %5.3f  %5.3f \n ',m, TimeSolver{dind}(index,1),TimeSolver{dind}(index,2));

        save Example3_5b_0707
    end

end
