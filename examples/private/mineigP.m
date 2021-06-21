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
    p0 = @(x,y) 1 - x.^2 - y.^2;
    p1 = @(x,y) x + x.*y - x.^3;
    p2 = @(x,y) 2.*x.^2.*y - x.*y - 2.*y.^3;
    for i = 1:row
        for j = 1:column
            r = input_r(i,j);
            theta = input_theta(i,j);
            x = r.*cos(theta);
            y = r.*sin(theta);
            P = ( p0(x,y).*eye(m) + p1(x,y).*A  + p2(x,y).*B ) .*r;
            lambda(i,j) = min(eig(P));
        end
    end
end

