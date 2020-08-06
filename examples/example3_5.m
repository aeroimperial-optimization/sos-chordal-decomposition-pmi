% =================================================================
% Code for Example 3.5 in the paper:
% Y. Zheng, G. Fantuzzi, Sum-of-squares chordal decomposition of polynomial 
% matrix inequalities
%
% NOTES: 
% 1. This code assumes that the option sos.csp in YALMIP has been modified,
%    or that one using AEROIMPERIAL_YALMIP -- see
%    https://github.com/aeroimperial-optimization/aeroimperial-yalmip
% 2. We plot the feasible set of a homogeneous polynomial matrix inequality
%    with chordal sparsity using a manual approach, detecting infeasibility and
%    exploiting symmetry. The code is ad-hoc but fast, based on YALMIP's
%    lmi/plot function.
% =================================================================

% Clean up
clear
close all

% Parameters -- EDIT HERE
d = 2;          % size of the matrix is 3*d
nTheta = 1e2;   % number of directions to optimize over  
                % 1e2 was used in the paper

% Set up P using an ugly loop
% Use a dirty trick to ensure that all variables x, y, z are coupled, so
% YALMIP ignores correlative sparsity in those variables when setting up
% the SOS programs with sparsity. This is needed because, since we use the
% scalarization u'*P(x,y,z)*u, not all terms [x*y, y*z, x*z]*u_i*u_j appear
% and YALMIP will exploit this!
sdpvar x y z
lambda = sdpvar(2,1);
sdpvar s
m = 3*d;
temp = s*(y^2*x^2+z^2*y^2+x^2*z^2);
a = [lambda(2)*x^4+y^4; lambda(2)*y^4+z^4; lambda(2)*z^4+x^4];
b = [x^2*y^2; y^2*z^2; z^2*x^2];
for i = 1:m-1
    if mod(i,2)
        c = lambda(1);
    else
        c = lambda(2);
    end
    ind = rem([i-1,i],3)+1;
    P(i,i) = a(ind(1));
    P(i+1,i+1) = a(ind(2));
    P(i,i+1) = b(ind(1))*c + temp; % dirty trick to force dense couplings in x,y,z
    P(i+1,i) = b(ind(1))*c + temp; % dirty trick to force dense couplings in x,y,z
end
clear c

% YALMIP options
opts = sdpsettings;
opts.cachesolvers = 1;
opts.sos.scale = 1;
opts.verbose = 0;

% Some useful variables
u = sdpvar(m,1);
theta = linspace(0,pi/2,nTheta);

% %% Dense SOS relaxation
% % Plot the feasible sets of the SOS approximations for all nu
% % Use scalarization. This code uses a manual trick to compile the SOS program 
% % once and change the objective function faster. This required some 
% % experimentation and the code looks obscure, but saves time
% nu = [1,2];     % exponent for (x^2 + y^2 + z^2)^nu weight
% runtime_csp0 = zeros(size(nu));
% opts.sos.csp = 0;
% feasset_csp0 = cell(length(nu),1);
% for j = 1:length(nu)
%     fprintf('Solving %i SOS programs (d=%i, nu=%i, CSP=%i): please wait...',nTheta,d,nu(j),opts.sos.csp)
%     % Compile the SOS constraints using YALMIP
%     % Create a low-level model like in lmi/plot
%     p = clean(u'*P*u * (x^2+y^2+z^2)^nu(j), 1e-12);
%     [F,h] = compilesos([sos(p),s==0],0,opts,[lambda;s]);
%     [model,recoverdata,diagnostic,internalmodel] = export(F,h,opts);
%     % Optimize objective calling MOSEK directly
%     feasset_csp0{j} = zeros(2,2*nTheta);
%     tstart = tic;
%     k = 0;
%     for i = 1:nTheta
%         model.prob.c(1:2) = [cos(theta(i)); sin(theta(i))];
%         [r, res] = mosekopt('minimize echo(0)',model.prob,model.param);
%         % Feasible or infeasible?
%         feasset_csp0{j}(:,i) = res.sol.itr.xx(1:2);
%         if strcmp(res.sol.itr.prosta, 'PRIMAL_AND_DUAL_FEASIBLE')
%             feasset_csp0{j}(:,i) = res.sol.itr.xx(1:2);
%             k = k + 1;
%         elseif any(strcmp(res.sol.itr.prosta, {'DUAL_INFEASIBLE', 'PRIMAL_INFEASIBLE'}))
%             feasset_csp0{j}(:,i) = [NaN; NaN];
%         end
%     end
%     runtime_csp0(j) = toc(tstart);
%     fprintf('done (%6.4f seconds)\n',runtime_csp0(j))
%     % Plot feasible set exploiting symmetry
%     % Remove NaNs, too...
%     feasset_csp0{j}(:,nTheta+1:end) = fliplr(feasset_csp0{j}(:,1:nTheta));
%     feasset_csp0{j}(1,nTheta+1:end) = -feasset_csp0{j}(1,nTheta+1:end);
%     feasset_csp0{j} = reshape(feasset_csp0{j}(~isnan(feasset_csp0{j})), 2, 2*nTheta);
% end


%% Sparse SOS relaxation
% Same as before but exploiting sparsity via the sos.csp option in YALMIP.
% This requires the dirty trick to get YALMIP to ignore correlative
% sparsity in the variables (x,y,z).
nu = [1];     % exponent for (x^2 + y^2 + z^2)^nu weight
runtime_csp1 = zeros(size(nu));
feasset_csp1 = cell(length(nu),1);
opts.sos.csp = 1;
for j = 1:length(nu)
    fprintf('Solving %i SOS programs (d=%i, nu=%i, CSP=%i): please wait...',nTheta,d,nu(j),opts.sos.csp)
    % Compile the SOS constraints using YALMIP
    % Create a low-level model like in lmi/plot
    p = clean(u'*P*u * (x^2+y^2+z^2)^nu(j), 1e-12);
    [F,h] = compilesos([sos(p),s==0],0,opts,[lambda;s]);
%     [F,h] = compilesos([sos(p)],0,opts,[lambda]);
    [model,recoverdata,diagnostic,internalmodel] = export(F,h,opts);
    model
    return
    % Optimize objective calling MOSEK directly
    feasset_csp1{j} = zeros(2,2*nTheta);
    tstart = tic;
    for i = 1:nTheta
        model.prob.c(1:2) = [cos(theta(i)); sin(theta(i))];
        [r, res] = mosekopt('minimize echo(0)',model.prob,model.param);
        % Feasible or infeasible?
        if strcmp(res.sol.itr.prosta, 'PRIMAL_AND_DUAL_FEASIBLE')
            feasset_csp1{j}(:,i) = res.sol.itr.xx(1:2);
        elseif any(strcmp(res.sol.itr.prosta, {'DUAL_INFEASIBLE', 'PRIMAL_INFEASIBLE'}))
            feasset_csp1{j}(:,i) = [NaN; NaN];
        end
    end
    runtime_csp1(j) = toc(tstart);
    fprintf('done (%6.4f seconds)\n',runtime_csp1(j))
    % Plot feasible set exploiting symmetry
    % Remove NaNs, too...
    feasset_csp1{j}(:,nTheta+1:end) = fliplr(feasset_csp1{j}(:,1:nTheta));
    feasset_csp1{j}(1,nTheta+1:end) = -feasset_csp1{j}(1,nTheta+1:end);
    feasset_csp1{j} = reshape(feasset_csp1{j}(~isnan(feasset_csp1{j})), 2, 2*nTheta);
end
