function plot_2D_field(X,Y,field,fieldtype,climmax,varnum,yshift,y_label)

field = field(:,:,varnum);

if nargin < 7
    yshift = 0;
end


if nargin < 8
    y_label = '$y$ [$\mu$m]';
end

vars = {'x','y','z'};
var = vars{varnum};

if strcmp(fieldtype,"E")
    c_label = sprintf('$|E_%s|^2$, [$(\\textrm{V}/\\textrm{m})^2$]',var);
elseif strcmp(fieldtype,"H")
    c_label = sprintf('$|H_%s|^2$, [$(\\textrm{A}/\\textrm{m})^2$]',var);
end

X = X*10^6;
Y = (Y- yshift)*10^6;

imagesc(X,Y,abs(field).^2)
set(gca,'YDir','normal')
axis equal
axis tight
if climmax > 0
    clim([0,climmax]);
end
hold on
c = colorbar;
c.Label.Interpreter = 'latex';
c.Label.String = c_label;
c.Label.FontSize = 12;
xlabel('$x$ [$\mu$m]','FontSize',14);
ylabel(y_label,'FontSize',14);
