
% =================================================================
% Draw Figure 3 in our paper 
% Y. Zheng, G. Fantuzzi, Sum-of-squares chordal decomposition of 
%                         polynomial matrix inequalities
% =================================================================

load Example3_5b
bnd_ctr = feasset_csp0;
bnd_dec = feasset_csp1;

% ============ plot feasible region ========================
set(0,'defaulttextinterpreter','latex')
set(0,'DefaultAxesColor','none')

% =========== no chordal decomposition =====================
figure;  
ColorBar = ['b','m','k','g'];
h = cell(2,1);
for k = 1:2
    patch(bnd_ctr{k}(1,:),bnd_ctr{k}(2,:),ColorBar(k),'FaceAlpha',0.08); hold on
    h{k} = plot(bnd_ctr{k}(1,:),bnd_ctr{k}(2,:),ColorBar(k),'linewidth',1.5);
end
h = legend([h{1},h{2}],'$$\nu = 1$$','$$\nu = 2$$','Location','Northwest');
set(h,'FontSize',12,'Interpreter','latex','box','off')


% Change axis
Xlim = [-1.5,1.5];
Ylim = [0,4.2];
xlim(Xlim);ylim(Ylim);
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

% font style
set(gca,'TickLabelInterpreter','latex','fontsize',12)
set(gca,'YTickLabel',{'0','1','2','3','4'},'YTick',[0 1 2 3 4]);
set(gca, 'Layer', 'top');

% labels
text(1.5,0.06,'$$\lambda_1$$','FontSize',14)
text(0.1,4.2,'$$\lambda_2$$','FontSize',14)
set(gcf,'Position',[100 100 300 300])
print(gcf,'Example3_5a','-painters','-dpng','-r600')


% =========== with chordal decomposition =====================
figure;  
ColorBar = ['b','m','g','k'];
h = cell(3,1);
for k = 2:3
    patch(bnd_dec{k}(1,:),bnd_dec{k}(2,:),ColorBar(k),'FaceAlpha',0.08); hold on
    h{k} = plot(bnd_dec{k}(1,:),bnd_dec{k}(2,:),ColorBar(k),'linewidth',1.5);
end
h = legend([h{2},h{3}],'$$\nu = 2$$','$$\nu = 3$$','Location','Northwest');
set(h,'FontSize',12,'Interpreter','latex','box','off')


% Change axis
Xlim = [-1.5,1.5];
Ylim = [0,4.2];
xlim(Xlim);ylim(Ylim);
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

% font style
set(gca,'TickLabelInterpreter','latex','fontsize',12)
set(gca,'YTickLabel',{'0','1','2','3','4'},'YTick',[0 1 2 3 4]);
set(gca, 'Layer', 'top');

% labels
text(1.5,0.06,'$$\lambda_1$$','FontSize',14)
text(0.1,4.2,'$$\lambda_2$$','FontSize',14)
set(gcf,'Position',[100 100 300 300])
print(gcf,'Example3_5b','-painters','-dpng','-r600')