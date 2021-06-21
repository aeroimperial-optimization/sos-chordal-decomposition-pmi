
% =================================================================
% Exact value of the integral function in Example 5-2 in the paper:
% Y. Zheng, G. Fantuzzi, Sum-of-squares chordal decomposition of polynomial 
% matrix inequalities
%
% =================================================================

function example5_2_integral
    load('data/example5_2_graphs.mat','A','B','Gsize')

    cost_integral = zeros(length(Gsize),1);
    for indm = 2%:length(Gsize)   % Different matris sizes
        m  = Gsize(indm);
        Am = A{indm};
        Bm = B{indm};
        p = @(r,theta) mineigP(m,Am,Bm,r,theta);
        tstart = tic;
        cost_integral(indm) = integral2(p,0,1,0,2*pi,'AbsTol',1e-7,'RelTol',1e-8); 
        tend = toc(tstart);
        fprintf('%d   %4.2f   %4.5f\n',m,tend,cost_integral(indm));
    end

    %save('example5_2_graphs.mat','G','A','B','Gsize','cost_integral')
end
