clear all
% Parameters -- EDIT HERE
nu = 3;    % exponent for (x^2 + y^2 + z^2)^nu weight
d = 2;     % size of the matrix is 3*d
plotset=0; % Plot feasible set? 0/1
% The code
sdpvar x y z
sdpvar k1 k2 s
k = [k1;k2];
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
    P(i,i+1) = b(ind(1))*c+temp; % dirty trick to force dense couplings in x,y,z
    P(i+1,i) = b(ind(1))*c+temp; % dirty trick to force dense couplings in x,y,z
end
clear c
sdisplay(replace(P,s,0))
u = sdpvar(m,1);
p = u'*P*u * (x^2+y^2+z^2)^nu;
opts = sdpsettings;
opts.cachesolvers = 1;
opts.sos.csp = 1;
% Plot feasible set
% Limit k2 by 4 otherwise unbounded
if plotset==1
    feasset = plot([sos(p),s==0,k2<=4],k,[],200);
    plot(feasset{1}(1,:),feasset{1}(2,:),'.-'); hold on
end
% Optimize
obj = k2-10*k1;
[sol,w,Q] = solvesos([sos(p),s==0],obj,opts,[k;s]);
kopt = value(k);
plot(kopt(1),kopt(2),'o')
