% =================================================================
% Plot inner approximation, Example 5-2 in the paper:
%
% Y. Zheng, G. Fantuzzi, Sum-of-squares chordal decomposition of polynomial 
% matrix inequalities
%
% =================================================================

% draw inner approximations
clear; close all
Dsparse = load('data/ex5_2_sparse.mat');  % Pre-computed results using sparse SOS formulation
% Ddense  = load('ex5_2_dense.mat');  % Pre-computed results using dense SOS formulation

% Clean up
yalmip clear

% Problem data  {x \in R^2, P(x) \geq 0}. 
sdpvar x y
p0 = @(x,y) 1 - x^2 - y^2;
p1 = @(x,y) x + x*y - x^3;
p2 = @(x,y) 2*x^2*y-x*y-2*y^3;

% Parameters
Deg = [2:6];            % half degree of SOS multipliers

% Parameters for figures 
N  = 100;
xg = linspace(-1,1,N); yg = linspace(-1,1,N);
ColorBar = lines(5); ColorBar = [0 0 0; ColorBar([1 2 4],:)];
Lwidth = 1;   % line width
Fontsize = 8;  % font size
FigSize = [8 3];

% set K = {1 - x1^2 - x2^2 \geq 0}; unit circle
th = linspace(0,2*pi,200);
xunit = 1 * cos(th) + 0;
yunit = 1 * sin(th) + 0;

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
    
    % Make figure of right size
    ff = figure('WindowStyle','normal');
    ff.Units = 'centimeters';
    ff.Position([3 4]) = FigSize;
    for indx = 1:length(Deg)
        subplot(1,length(Deg),indx)
        
        % Plot the unit circle 
        plot(xunit,yunit,'k:','linewidth',Lwidth);hold on
        
        % plot the real boundary 
        [bnd0,h0]=contour(xg_mesh,yg_mesh,Eig_grid,[0 0],'color',ColorBar(1,:),'linewidth',Lwidth); 

        % Dense SOS
%         monobasis = monolist([x,y],2*Deg(indx));
%         if ~isempty(gStandard{indx,indm})  % has a solution
%         p = gStandard{indx,indm}'*monobasis;  % { x\in R^2, p(x) >=0}
%         pgrid = zeros(N,N);
%         for i = 1:N
%             for j = 1:N
%                 pgrid(i,j) = replace(p,[x y],[xg(i) yg(j)]);
%                 if xg(i)^2 +  yg(j)^2 > 1
%                     pgrid(i,j) = -100;  % outside the unit circle
%                 end
%             end
%         end
%         hold on; [bnd1,h1] = contour(xg_mesh,yg_mesh,pgrid',[0 0],ColorBar(2),'linewidth',Lwidth);
%         end

        % Sparse SOS
        pgrid = zeros(N,N);
        powers = Dsparse.exponents{indx,indm};
        for i = 1:size(powers,1)
            pgrid = pgrid + Dsparse.gSparse{indx,indm}(i) .* xg_mesh.^powers(i,1) .* yg_mesh.^powers(i,2);
        end
        pgrid(xg_mesh.^2+yg_mesh.^2>1) = -100;
        hold on; [bnd2,h2] = contour(xg_mesh,yg_mesh,pgrid,[0 0], 'color',ColorBar(3,:),'linewidth',Lwidth);
        
        % Format
        axis square
        axis([-1 1 -1 1])
        xlabel('$x_1$','Interpreter','latex','FontSize',Fontsize)
        if indx==1; ylabel('$x_2$','Interpreter','latex','FontSize',Fontsize); end
        set(gca,'TickLabelInterpreter','latex','fontsize',Fontsize);
    end

    % Print
    ff.Position([3 4]) = FigSize;
    fname = ['inner_approx',num2str(Dsparse.Gsize(indm))];
%     print(gcf,fname,'-painters','-depsc','-r600')
end


