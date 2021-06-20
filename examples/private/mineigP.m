function lambda = mineigP(m,A,B,input_r,input_theta)
% minimum eigenvalue
%
% p0 = 1 - x^2 - y^2;
% p1 = x + x*y - x^3;
% p2 = 2*x^2*y-x*y-2*y^3;
%
% P = p0*I + p1*A + p2*B
% min eig(P)

    [row,column]  = size(input_r);
    lambda = zeros(size(input_r));
    for i = 1:row
        for j = 1:column
            r = input_r(i,j);
            theta = input_theta(i,j);
            P = (1 - (r*sin(theta))^2 - (r*cos(theta))^2)*eye(m) ...
                + (r*sin(theta) + r*sin(theta)*r*cos(theta) - (r*sin(theta))^3)*A ...
                + (2*(r*sin(theta))^2*r*cos(theta) - r*sin(theta)*r*cos(theta)-2*(r*cos(theta))^3)*B;
            lambda(i,j) = min(eig(P));
        end
    end
end

