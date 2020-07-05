
% =================================================================
% Code for Example 3-5 in our paper
% Y. Zheng, G. Fantuzzi, Sum-of-squares chordal decomposition of 
%                         polynomial matrix inequalities
%
% Assume the option sos.csp in Yalmip has been modified
% =================================================================

clc;clear;
sdpvar x y z        % sos variables
v = sdpvar(3,1);    % scalarization  
sdpvar a b          % parameters

P = [x^4+y^4,    -a*x^2*y^2,   0; 
     -a*x^2*y^2, a*y^4+b*z^4, -b*y^2*z^2; 
     0,          -b*y^2*z^2,    x^4+z^4];
opts =  sdpsettings('solver','mosek');

%% standard SOS
opts.sos.csp = 0;
bnd_ctr  = cell(3,1);
for k = 1:4
    F = sos((x^2+y^2+z^2)^(k-1)*v'*P*v);               % The current Yalmip only exploit csp in the scalar case
    tmp = plot(F,[a,b],'m',100,opts);                  %% boundary points
    bnd_ctr{k} = tmp{1};
end

%% chordal decomposition
opts.sos.csp = 1;
bnd_dec  = cell(3,1);
for k = 1:4
    F = sos((x^2+y^2+z^2)^(k-1)*v'*P*v);     % The current Yalmip only exploit csp in the scalar case
    tmp = plot(F,[a,b],'m',100,opts);        %% boundary points
    bnd_dec{k} = tmp{1};
end

%save Example3_5


