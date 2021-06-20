% =================================================================
% Code for Example 5-2 in the paper:
% Y. Zheng, G. Fantuzzi, Sum-of-squares chordal decomposition of polynomial 
% matrix inequalities
%
% Assume the option sos.csp in Yalmip has been modified
% =================================================================

function example5_2
    % Clean up
    clear
    yalmip clear

    % load the graph data
    load('data/example5_2_graphs.mat','G','A','B','Gsize')  % pre-generated chordal graphs and data

    % Problem data  {x \in R^2, P(x) \geq 0}. 
    sdpvar x y
    p0 = 1 - x^2 - y^2;
    p1 = x + x*y - x^3;
    p2 = 2*x^2*y-x*y-2*y^3;

    % Parameters
    Deg = [2,3,4];                    % half degree of SOS multipliers
    
    % initializing the results
    gStandard    = cell(length(Deg),length(Gsize));
    objStandard  = zeros(length(Deg),length(Gsize));
    timeStandard = zeros(length(Deg),length(Gsize));

    gSparse    = cell(length(Deg),length(Gsize));
    objSparse  = zeros(length(Deg),length(Gsize));
    timeSparse = zeros(length(Deg),length(Gsize));
    cons_time  = zeros(length(Deg),length(Gsize));
    
    %load Example5_5.mat
    for indm = 1:length(Gsize)
        m         = Gsize(indm);
        fprintf('Graph size: %d \n\n',m);
        
        P = p0*eye(m)+ p1*A{indm} + p2*B{indm};

        % region B that contains P, which is a unit circle in our case
        g0 = 1 - x^2 - y^2;

        v   = sdpvar(m,1);        % scalarization
        % --------------------------------------------------
        % Standard SOS optimization
        % --------------------------------------------------
        memFlag = 0;  % out of memory for Standard SOS?
        for indx = 1:length(Deg)  % long over degree

            % ----------------------------------------------
            % set up the cost function
            % ----------------------------------------------

            % objective function g(x) >= 0 
            monobasis = monolist([x,y],2*Deg(indx));
            coeff     = sdpvar(length(monobasis),1);
            g         = coeff'*monobasis;

            % Calculate the moment sequence -- integral part  int xy dxdy
            const_integral = zeros(length(coeff),1);
            for i = 1:length(coeff)   
                tmp = degree(monobasis(i),[x y]);
                p = @(r,theta) (r.*sin(theta)).^(tmp(1)).*(r.*cos(theta)).^(tmp(2)); % polar form
                const_integral(i) = integral2(p,0,1,0,2*pi);  % polar form over a unit circle
            end
            cost = const_integral'*coeff;

            % ----------------------------------------------
            % set up the SOS constraint
            % ----------------------------------------------
            monobasis = monolist([x,y],Deg(indx));
            monobasis_Matrix = kron(eye(m),monobasis);
            N  = length(monobasis);
            Q0 = sdpvar(m*N); S0 = monobasis_Matrix'*Q0*monobasis_Matrix;   % SOS matrix multipler
            s0 = v'*S0*v;

            monobasis = monolist([x,y],Deg(indx) - 1);
            monobasis_Matrix = kron(eye(m),monobasis);
            N  = length(monobasis);
            Q1 = sdpvar(m*N); S1 = monobasis_Matrix'*Q1*monobasis_Matrix;   % SOS matrix multipler
            s1 = v'*S1*v;

            F = [sos(s0),sos(s1)];
            F = [F, sos(v'*(P-g*eye(m))*v - s0 - s1*g0)];

            % ----------------------------------------------
            % Call a solver 
            % ----------------------------------------------
            sosQvar = [vec(Q0);vec(Q1)];
            opts = sdpsettings('sos.csp',0,'solver','mosek');  % no correlative sparsity
            if m >= 30 && Deg(indx) == 4  % these choices --> out of memory on our computer
                memFlag = 1;
            end
            if m >= 35 && Deg(indx) >= 3 % these choices --> out of memory on our computer
                memFlag = 1;
            end
            if memFlag == 0
            try
                sol =  solvesos(F,-cost,opts,[sosQvar;coeff]);
                % record solution
                    gStandard{indx,indm}    = value(coeff);
                    objStandard(indx,indm)  = value(cost);
                    timeStandard(indx,indm) = sol.solvertime;
            catch
                memFlag = 1;   % don't need to try higher degrees
                warning('out of memory')
            end
            end

        end

        % --------------------------------------------------
        % Sparse SOS optimization
        % --------------------------------------------------

        for indx = 1:length(Deg)  % long over degree

            % ----------------------------------------------
            % set up the cost function
            % ----------------------------------------------

            % objective function g(x) >= 0 
            monobasis = monolist([x,y],2*Deg(indx));
            coeff     = sdpvar(length(monobasis),1);
            g         = coeff'*monobasis;

            % Calculate the moment sequence -- integral part  int xy dB
            const_integral = zeros(length(coeff),1);
            for i = 1:length(coeff)   
                tmp = degree(monobasis(i),[x y]);
                p = @(r,theta) (r.*sin(theta)).^(tmp(1)).*(r.*cos(theta)).^(tmp(2)); % polar form
                const_integral(i) = integral2(p,0,1,0,2*pi);  % polar form over a unit circle
            end
            cost = const_integral'*coeff;

            % ----------------------------------------------
            % set up the SOS constraint
            % ----------------------------------------------
              clear s0 s1 Q0 Q1
              clique = cliquesFromSpMatD(G{indm});
              tic    % the following construction seems very time consuming
              s0 = cell(clique.NoC,1); Q0 = cell(clique.NoC,1);
              s1 = cell(clique.NoC,1); Q1 = cell(clique.NoC,1);

              F = [];           % SOS constraint
              temp = 0;
              sosQvar = [];     % Gram matrix Q
              for t = 1:clique.NoC  
                  % s0
                  monobasis = monolist([x,y],Deg(indx));
                  monobasis_Matrix = kron(eye(clique.NoElem(t)),monobasis);
                  N = length(monobasis);
                  Q0{t} = sdpvar(clique.NoElem(t)*N);
                  S = monobasis_Matrix'*Q0{t}*monobasis_Matrix;   % SOS matrix
                  s0{t} = v(clique.Set{t},1)'*S*v(clique.Set{t},1);

                  % s1
                  monobasis = monolist([x,y],Deg(indx)-1);
                  monobasis_Matrix = kron(eye(clique.NoElem(t)),monobasis);
                  N = length(monobasis);
                  Q1{t} = sdpvar(clique.NoElem(t)*N);
                  S = monobasis_Matrix'*Q1{t}*monobasis_Matrix;   % SOS matrix
                  s1{t} = v(clique.Set{t},1)'*S*v(clique.Set{t},1);

                  F = [F,sos(s0{t}),sos(s1{t})];
                  temp = temp + s0{t} + s1{t}*g0;
                  sosQvar = [sosQvar;vec(Q0{t});vec(Q1{t})];
              end
              cons_time(indx,indm) = toc;

              F = [F, sos(v'*(P-g*eye(m))*v - temp)];

            % ----------------------------------------------
            % Call a solver 
            % ----------------------------------------------
            opts = sdpsettings('sos.csp',1,'solver','mosek');  % exploit correlative sparsity
            sol =  solvesos(F,-cost,opts,[sosQvar;coeff]);

            % record solution
            gSparse{indx,indm}    = value(coeff);
            objSparse(indx,indm)  = value(cost);
            timeSparse(indx,indm) = sol.solvertime;
        end
        %save('example5_2_results.mat','Gsize','G','A','B','gStandard','objStandard','timeStandard','gSparse','objSparse','timeSparse','cons_time');
    end
end