% =================================================================
% Generating random chordal patterns, Example 5-2 in the paper:
% Y. Zheng, G. Fantuzzi, Sum-of-squares chordal decomposition of polynomial 
% matrix inequalities
%
% =================================================================

Gsize = 5:2:10;  % 10:5:40 was used in the paper
maxClique = 5;   % bound on the largest maximal cliques
G = [];
A = [];
B = [];

for indm = 1:length(Gsize)
    m         = Gsize(indm);
    G{indm}   = chordalGen(m,maxClique);         % Chordal graph with largest clique size <= maxClique
    tmp = rand(m);
    tmp = (tmp + tmp')/2;                        % symmetric matrix 
    A{indm} = (spones(G{indm}) - eye(m)).*tmp;   % zeros on its diagonal

    tmp = rand(m);
    tmp = (tmp + tmp')/2;      % symmetric
    B{indm} = (spones(G{indm}) - eye(m)).*tmp;   % zeros on its diagonal
end

%save('data/example5_2_graphs_new.mat','G','A','B','Gsize');