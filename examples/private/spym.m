function spym(X,varargin)

% SPYM.m spy matrix, formatted
%
% 

% Spy

if isa(X,'sdpvar')
    X = any(X);
end
[n,m] = size(X);
spy(sparse(X),varargin{:});
grid on

% Format
h = gca;
axis([0.5 m+0.5 0.5 n+0.5])
h.YTick = 0.5:1:n+0.5;
h.XTick = 0.5:1:m+0.5;
h.YTickLabel = {};
h.XTickLabel = {};
h.XLabel.String = '';
h.GridLineStyle = '-';
h.TickLength = [0 0];
h.Children.MarkerFaceColor = h.Children.Color;



