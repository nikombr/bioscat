function segment = setup_nanostructures(segx,segy,alpha)

d = segx(2)-segx(1);
r = min((max(segx) - min(segx)),(max(segy) - min(segy)))/2;


%segy = ones(size(segy))*max(segy);

x_top    = segx;
y_top    = segy;
ys = 0:d:segy(end);
ys = linspace(0,segy(end),length(ys));
x_right  = 0*ys+max(segx);
y_right  = flip(ys);
x_bottom = flip(segx);
y_bottom = 0*segy;
ys = 0:d:segy(1);
ys = linspace(0,segy(1),length(ys));
x_left = 0*ys+min(segx);
y_left = ys;

x = {x_top, x_right, x_bottom, x_left};
y = {y_top, y_right, y_bottom, y_left};

%ys = 0:d:segy(end);
%x = [segx 0*ys+max(segx)];
%y = [segy flip(ys)];
%x = [x flip(segx)];
%y = [y 0*segy];
%ys = 0:d:segy(1);
%x = [x 0*ys+min(segx)];
%y = [y ys];

X = {x_top, x_right, x_bottom, x_left};
Y = {y_top, y_right, y_bottom, y_left};

% Compute exterior points
x = [x_top x_right x_bottom x_left];
y = [y_top y_right y_bottom y_left];

x1  = [x(end) x(1:end-1)];
y1  = [y(end) y(1:end-1)];
x2 = [x(2:end) x(1)];
y2 = [y(2:end) y(1)];
x_diff = x1 - x2;
y_diff = y1 - y2;
diff = [x_diff; y_diff];
diff = diff./vecnorm(diff)*alpha;

x_ext = x + diff(2,:);
y_ext = y - diff(1,:);

n_x = x_ext(1:length(x_top));
n_y = y_ext(1:length(y_top));
n_x = n_x(2:end-1);
n_y = n_y(2:end-1);

% Compute interior points
x_int = [];
y_int = [];
for k = 1:4
    xseg = X{k};
    yseg = Y{k};
    x = xseg(3:end-2);
    y = yseg(3:end-2);
    
    x1  = xseg(2:end-3);
    y1  = yseg(2:end-3);
    x2  = xseg(4:end-1);
    y2  = yseg(4:end-1);
    x_diff = x1 - x2;
    y_diff = y1 - y2;
    diff = [x_diff; y_diff];
    diff = diff./vecnorm(diff)*alpha;

    xval = x - diff(2,:);
    yval = y + diff(1,:);

    x_int = [x_int xval];
    y_int = [y_int yval];

end

% Remove first and last point in each vector to avoid corners
x_top    = x_top(2:end-1);
y_top    = y_top(2:end-1);
x_right  = x_right(2:end-1);
y_right  = y_right(2:end-1);
x_bottom = x_bottom(2:end-1);
y_bottom = y_bottom(2:end-1);
x_left   = x_left(2:end-1);
y_left   = y_left(2:end-1);

% Save data in struct
segment = struct;
segment.x_top    = x_top';
segment.y_top    = y_top';
segment.n_x    = n_x';
segment.n_y    = n_y';
segment.x_right  = x_right';
segment.y_right  = y_right';
segment.x_bottom = x_bottom';
segment.y_bottom = y_bottom';
segment.x_left = x_left';
segment.y_left = y_left';
segment.x = {segment.x_top, segment.x_right, segment.x_bottom, segment.x_left};
segment.y = {segment.y_top, segment.y_right, segment.y_bottom, segment.y_left};
segment.x_int = x_int;
segment.y_int = y_int;
segment.x_ext = x_ext;
segment.y_ext = y_ext;