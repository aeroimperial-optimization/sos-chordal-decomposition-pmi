function [x, fval] = example3_6_fmincon(m,plotAll)
% Exact answer for example 3.6 using fmincon
% Inputs:
% * m: size of the matrix P
% * plotAll: true to make a contour plot of the minmum eigenvalue of P on
%            the semialgebraic set K. The code to do this is really
%            inefficient, but hey...it's simple

% Plot?
if nargin<2
    plotAll = 0;
end

% Make contour plot?
if plotAll
    [xx,yy] = meshgrid(-1:0.01:1,-1:0.01:1);
    E = NaN(size(xx));
    for j =1:size(xx,2)
        for i = 1:size(xx,1)
            if (xx(i,j)^2-1 <= 0) && (yy(i,j)^2 - xx(i,j)^2 <= 0)
                E(i,j) = OBJ([xx(i,j),yy(i,j)], m);
            end
        end
    end
    contourf(xx,yy,E);
    xlabel('x')
    ylabel('y')
    colorbar
end

% Try fmincon (carefully selected initial condition)
x0 = [0.9, -0.5];
options = optimoptions('fmincon');
options.FunctionTolerance = 1e-14;
options.StepTolerance = 1e-14;
options.OptimalityTolerance = 1e-14;
[x,fval] = fmincon(@(x)OBJ(x,m),x0,[],[],[],[],[],[],@(x)NONLCON(x),options);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E = OBJ(x,m)
% Minimum eigenvalue of P
a = 10 + x(2).^3 - x(1).^4;            % diagonal entry
b = x(1) + x(1).*x(2) - x(1).^3;       % off-diagonal entry

% Evaluate matrix
P = diag( repmat(a,m,1) );
P(2:end,1) = b;
P(1,2:end) = b;

% Min eig
E = min(eig(P));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C,EQ] = NONLCON(x)
% constraint x to the semialgebraic set. Exploit symmetry so consider only
% the positive part of the bow-tie set.
EQ = [];
C = [x(1)^2-1; x(2)^2-x(1)^2];
end

