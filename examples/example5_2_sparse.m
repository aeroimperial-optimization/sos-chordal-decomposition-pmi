% =================================================================% Code for Example 5-2 in the paper:
% Y. Zheng, G. Fantuzzi, Sum-of-squares chordal decomposition of polynomial
% matrix inequalities
%
% This script implements the sparse SOS relaxations, using the scalarized
% version of the decomposition results in the paper.
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
Deg = 10;

% Initalize containers for the results
gSparse    = cell(length(Deg),length(Gsize));
exponents    = cell(length(Deg),length(Gsize));
objSparse  = zeros(length(Deg),length(Gsize));
timeSparse = zeros(length(Deg),length(Gsize));

for indm = 2%length(Gsize)
    m = Gsize(indm);
    clique = cliquesFromSpMatD(G{indm}); % cliques of graph
    fprintf('     Graph size: %d \n',m)
    fprintf('Max clique size: %d \n\n',max(cellfun(@length, clique.Set)))
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
        exponents{indx,indm} = getexponentbase(g,[x,y]);
        exponents{indx,indm} = full(exponents{indx,indm});
        moment = zeros(size(exponents{indx,indm},1),1);
        for i = 1:size(exponents{indx,indm},1)
            % Integrate x^a*y^b over unit disk in polar coordinates
            tmp = exponents{indx,indm}(i,:);
            p = @(r,theta) (r.*cos(theta)).^(tmp(1)) .* (r.*sin(theta)).^(tmp(2)) .*r;
            moment(i) = integral2(p,0,1,0,2*pi);
        end
        cost = moment.'*gc;
        
        % Set up the SOS constraints, exploiting sparsity
        CNSTR = [];
        poly = v'*( P - g*eye(m) )*v;
        symmetries = cell(clique.NoC, 1);
        for k = 1:clique.NoC
            [s{k}, sc{k}] = multidegpoly({[x y], v(clique.Set{k})}, [2*Deg(indx)-2, 2], [0,2]);
            CNSTR = CNSTR + sos(s{k});
            poly = poly - g0 * s{k};
            
            % The constraints have symmetries! It shouldn't make a
            % difference, but we specify them anyway and let YALMIP figure
            % out if they are useful or not. We specify them because, when
            % m>14, YALMIP turns off automatic symmetry detection
            symmetries{k} = [0; 0; ones(length(clique.Set{k}), 1)];
        end
        CNSTR = [CNSTR; sos(poly)];
        symmetries{end+1} = [0; 0; ones(m, 1)];  
        
        % Solve, exploiting sparsity
        opts = sdpsettings('sos.csp',1);
        params = [vertcat(sc{:}); gc];
        sol = solvesos(CNSTR,-cost,opts,params,[],symmetries);
        gSparse{indx,indm}    = value(gc);
        objSparse(indx,indm)  = value(cost);
        timeSparse(indx,indm) = sol.solvertime;
    end
%     save(sprintf('./data/ex5_2_sparse_graph%i.mat',indm),'Gsize','G','A','B','gSparse','objSparse','timeSparse','exponents');
end
