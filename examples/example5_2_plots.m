% =================================================================
% Plot inner approximation (Figure 6), Example 5-2 in the paper:
%
% Y. Zheng, G. Fantuzzi, Sum-of-squares chordal decomposition of polynomial
% matrix inequalities
%
% =================================================================

% draw inner approximations
clear; close all
Dsparse = load('data/ex5_2_sparse_new.mat');  % Pre-computed results using sparse SOS formulation
Ddense  = load('data/ex5_2_dense_new.mat');  % Pre-computed results using dense SOS formulation

% Clean up
yalmip clear

% Problem data  {x \in R^2, P(x) \geq 0}.
p0 = @(x,y) 1 - x^2 - y^2;
p1 = @(x,y) x + x*y - x^3;
p2 = @(x,y) 2*x^2*y-x*y-2*y^3;

% Parameters
Deg = [2,3,4];            % half degree of SOS multipliers

% Parameters for figures
N  = 100;
xg = linspace(-1,1,N); yg = linspace(-1,1,N);
ColorBar = lines(5); ColorBar = [0 0 0; ColorBar([1 2 4],:)];
Lwidth = 1;   % line width
Fontsize = 8;  % font size
FigSize = [15 8];

% set K = {1 - x1^2 - x2^2 \geq 0}; unit circle
th = linspace(0,2*pi,200);
xunit = 1 * cos(th) + 0;
yunit = 1 * sin(th) + 0;

% Make figure of right size
ff = figure('WindowStyle','normal');
ff.Units = 'centimeters';
ff.Position([3 4]) = FigSize;
pltInd = reshape(1:18,6,3)';
shift = 0;
for indm = 2:length(Dsparse.Gsize)   % Different matris sizes
    
    % ------------------------------------------------------------------------
    % Figures: Plot the contour = 0
    % ------------------------------------------------------------------------
    % Plot the real PSD region  {x \in R^2, P(x) \geq 0}
    [xg_mesh,yg_mesh] = meshgrid(xg,yg);
    Eig_grid = zeros(N,N);
    for i = 1:numel(xg_mesh)
        x = xg_mesh(i);
        y = yg_mesh(i);
        P = p0(x,y).*eye(Dsparse.Gsize(indm)) + p1(x,y)*Dsparse.A{indm} + p2(x,y)*Dsparse.B{indm};
        Eig_grid(i) = min(eig(P));
    end
    
    for indx = 1:length(Deg)
        shift = shift+1;
        
        % Exact stuff
        subplot(length(Deg),6,pltInd(shift))
        plot(xunit,yunit,'k:','linewidth',Lwidth); hold on
        [bnd0,h0]=contour(xg_mesh,yg_mesh,Eig_grid,[0 0],'color',ColorBar(1,:),'linewidth',Lwidth);
        
        % Dense SOS
        pgrid = zeros(N,N);
        powers = Ddense.exponents{indx,indm};
        if ~isempty(Ddense.gStandard{indx,indm})  % value computed
            for i = 1:size(powers,1)
                pgrid = pgrid + Ddense.gStandard{indx,indm}(i) .* xg_mesh.^powers(i,1) .* yg_mesh.^powers(i,2);
            end
            pgrid(xg_mesh.^2+yg_mesh.^2>1) = NaN;
            hold on; [bnd1,h1] = contour(xg_mesh,yg_mesh,pgrid,[0 0], 'color',ColorBar(2,:),'linewidth',Lwidth);
            Z = pgrid; Z(Z<0) = NaN; Z(Z>=0) = 1;
            colormap(ColorBar(2,:))
            p = pcolor(xg_mesh,yg_mesh,Z);
            p.LineStyle = 'none';
            p.FaceAlpha = 0.25;
        end
        
        
        % Sparse SOS
        pgrid = zeros(N,N);
        powers = Dsparse.exponents{indx,indm};
        for i = 1:size(powers,1)
            pgrid = pgrid + Dsparse.gSparse{indx,indm}(i) .* xg_mesh.^powers(i,1) .* yg_mesh.^powers(i,2);
        end
        pgrid(xg_mesh.^2+yg_mesh.^2>1) = -100;
        hold on; [bnd2,h2] = contour(xg_mesh,yg_mesh,pgrid,[0 0], 'color',ColorBar(3,:),'linewidth',Lwidth);
        h2.LineStyle = '-';
        
        % Format
        axis square
        axis([-1 1 -1 1])
        ax = gca;
        ax.XTickLabels = [];
        set(gca,'TickLabelInterpreter','latex','fontsize',Fontsize);
        if indx==length(Deg)
            xlabel('$x_1$','Interpreter','latex','FontSize',Fontsize);
            ax.XTickLabels = [-1 0 1];
        end
        if indm==2; yPos(indx)=ax.Position(2); end
        if indm==2
            ylabel('$x_2$','Interpreter','latex','FontSize',Fontsize);
        elseif indm>2
            ax.YTickLabel = [];
        end
        t = title(sprintf('$m=%i$, $d=%i$',Dsparse.Gsize(indm),Deg(indx)));
        t.Interpreter = 'latex';
        t.FontSize = Fontsize-2;
        ax.Position([3 4]) = 1.1*ax.Position([3 4]);
        ax.Position(2) = yPos(indx);
    end
end

% Print
drawnow
pause(0.5)
ff.Position([3 4]) = FigSize;
ff.PaperUnits = 'centimeters';
ff.PaperSize = FigSize;
ff.PaperPosition = [0 0 FigSize];
fname = ['gio_inner_approx'];
%print(gcf,fname,'-painters','-depsc','-r300')


