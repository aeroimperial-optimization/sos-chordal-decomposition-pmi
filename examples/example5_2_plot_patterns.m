% =================================================================
% Plot chordal patterns (Figure 5), Example 5-2 in the paper:
% Y. Zheng, G. Fantuzzi, Sum-of-squares chordal decomposition of polynomial 
% matrix inequalities
%
% =================================================================

close all
load('data/example5_2_graphs.mat','G','Gsize');  % Pre-generated chordal graphs
for indg = 2:7
    figure;
    spym(G{indg},'ks',5);               
    set(gcf,'Position',[100 100 300 300])
    fname = ['chordal',num2str(Gsize(indg))];
    %print(gcf,fname,'-painters','-depsc','-r600')
end

