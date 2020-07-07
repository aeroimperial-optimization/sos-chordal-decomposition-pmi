clc; clear all; 

% Parameters -- EDIT HERE

%% Generating the polynomial matrix

sdpvar x y z
sdpvar k1 k2 s

k = [k1;k2];
d = 2;     % size of the matrix is 3*d
m = 3*d;
temp = s*(y^2*x^2+z^2*y^2+x^2*z^2);
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

%sdisplay(replace(P,s,0))

sdisplay(P)


%% Plot feasible set
% Limit k2 by 4 otherwise unbounded
v = sdpvar(m,1);
opts = sdpsettings;
opts.cachesolvers = 1;


%% standard SOS
opts.sos.csp = 0;
bnd_ctr  = cell(3,1);
nu = [2,3,4];
for ind = 1:3
    F = sos((x^2+y^2+z^2)^(nu(ind)-1)*v'*P*v);               % The current Yalmip only exploit csp in the scalar case
    tmp = plot([F,k2<=4],[k1,k2],'m',100,opts);                  %% boundary points
    bnd_ctr{ind} = tmp{1};
end


ColorBar = ['b','m','k','g'];
figure;
for k = 1:3
    patch(bnd_ctr{k}(1,:),bnd_ctr{k}(2,:),ColorBar(k),'FaceAlpha',0.08); hold on
    h{k} = plot(bnd_ctr{k}(1,:),bnd_ctr{k}(2,:),ColorBar(k),'linewidth',1.5);
end

%% chordal decomposition
opts.sos.csp = 1;
bnd_dec  = cell(3,1);
for ind = 1:3
    F = sos((x^2+y^2+z^2)^(nu(ind)-1)*v'*P*v);     % The current Yalmip only exploit csp in the scalar case
    tmp = plot([F,k2<=4],[k1,k2],'m',100,opts);        %% boundary points
    bnd_dec{ind} = tmp{1};
end


figure;
for k = 1:3
    patch(bnd_dec{k}(1,:),bnd_dec{k}(2,:),ColorBar(k),'FaceAlpha',0.08); hold on
    h{k} = plot(bnd_dec{k}(1,:),bnd_dec{k}(2,:),ColorBar(k),'linewidth',1.5);
end

save Example3_5_v2




