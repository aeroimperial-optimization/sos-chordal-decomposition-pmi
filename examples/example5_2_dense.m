% =================================================================
% Code for Example 5.2 in the paper:
% Y. Zheng, G. Fantuzzi, "Sum-of-squares chordal decomposition of polynomial
% matrix inequalities", 2021
%
% This script implements the standard SOS relaxations, using the scalarized
% formulation.
%
% NOTE: To run this script, please install the following fork of YALMIP:
% https://github.com/aeroimperial-optimization/aeroimperial-yalmip
% =================================================================

% Clean up
clear, clc
yalmip clear

% Load pre-generated chordal graphs and data
load('data/example5_2_graphs.mat','G','A','B','Gsize')

% Parameters: half-degree of SOS multipliers
Deg = [2,3,4];

% Initalize containers for the results
gStandard    = cell(length(Deg),length(Gsize));
objStandard  = zeros(length(Deg),length(Gsize));
timeStandard = zeros(length(Deg),length(Gsize));

for indm = 1%:length(Gsize)
    m         = Gsize(indm);
    fprintf('Graph size: %d \n\n',m);
    memFlag = 0;  % out of memory for Standard SOS?
    for indx = 1:length(Deg)  % long over degree
        
        % Clear yalmip variables (avoid build-up of symbolic variables that
        % slows down yalmip
        yalmip clear
        sdpvar x y
        p0 = 1 - x^2 - y^2;
        p1 = x + x*y - x^3;
        p2 = 2*x^2*y - x*y - 2*y^3;
        P = p0*eye(m)+ p1*A{indm} + p2*B{indm};
        g0 = 1 - x^2 - y^2;
        v   = sdpvar(m,1);
        
        % Set up the cost function
        % Must integrate all moments x^a*y^b over the unit disk; work in
        % polar coordinates for simplicity, and use numerical integration
        % even though analytical expressions exist!
        % Note: area element in polar coordinates is r*dr*d\theta
        [g,gc] = polynomial([x,y],2*Deg(indx));
        exponents = getexponentbase(g,[x,y]);
        exponents = full(exponents);
        moment = zeros(size(exponents,1),1);
        for i = 1:size(exponents,1)
            % Integrate x^a*y^b over unit disk in polar coordinates
            tmp = exponents(i,:);
            p = @(r,theta) (r.*cos(theta)).^(tmp(1)) .* (r.*sin(theta)).^(tmp(2)) .*r;
            moment(i) = integral2(p,0,1,0,2*pi);
        end
        cost = moment.'*gc;
        
        % Set up the SOS constraints
        [s1, s1c] = multidegpoly({[x y], v}, [2*Deg(indx)-2, 2], [0,2]);
        s0 = v'*( P - g*eye(m) )*v - s1*g0; % The actual SOS constraint
        F = [sos(s0), sos(s1)];
        
        % Solve without sparsity exploitation
        opts = sdpsettings('sos.csp',0);
        if m >= 30 && Deg(indx) == 4  % these choices --> out of memory on our computer
            memFlag = 1;
        end
        if m >= 35 && Deg(indx) >= 3 % these choices --> out of memory on our computer
            memFlag = 1;
        end
        if memFlag == 0
            try
                sol =  solvesos(F,-cost,opts,[s1c; gc]);
                gStandard{indx,indm}    = value(gc);
                objStandard(indx,indm)  = value(cost);
                timeStandard(indx,indm) = sol.solvertime;
            catch
                memFlag = 1;   % don't need to try higher degrees
                warning('out of memory')
            end
        end
        
    end
end