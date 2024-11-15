function segment = setup_nanostructures(segx,segy,alpha)

d = segx(2)-segx(1);
r = min((max(segx) - min(segx)),(max(segy) - min(segy)))/2;
minNumSteps = 10;
numstepstart = max(minNumSteps,ceil(segy(1)/d));
numstepend = max(minNumSteps,ceil(segy(end)/d));

x_top    = segx;
y_top    = segy;
ys = linspace(0, segy(end), numstepend+1);
x_right  = 0*ys+max(segx);
y_right  = flip(ys);
x_bottom = flip(segx);
y_bottom = 0*segy;
ys = linspace(0, segy(1), numstepstart+1);
x_left = 0*ys+min(segx);
y_left = ys;

X = {x_top, x_right, x_bottom, x_left};
Y = {y_top, y_right, y_bottom, y_left};

% Compute exterior points
x = [x_top(1:end-1) x_right(1:end-1) x_bottom(1:end-1) x_left(1:end-1)];
y = [y_top(1:end-1) y_right(1:end-1) y_bottom(1:end-1) y_left(1:end-1)];

x1  = [x(end) x(1:end-1)];
y1  = [y(end) y(1:end-1)];
x2 = [x(2:end) x(1)];
y2 = [y(2:end) y(1)];
x_diff = x1 - x2;
y_diff = y1 - y2;
diff = [x_diff; y_diff];
diff = diff./vecnorm(diff);

x_ext = x + diff(2,:)*alpha;
y_ext = y - diff(1,:)*alpha;

n_x = diff(2,1:length(x_top));
n_y = -diff(1,1:length(y_top));
n_x = n_x(2:end-1);
n_y = n_y(2:end-1);

% Compute interior points
x_int = [];
y_int = [];
for k = 1:4
    xseg = X{k};
    yseg = Y{k};
    %x = xseg(3:end-2);
    %y = yseg(3:end-2);
    %x1  = xseg(2:end-3);
    %y1  = yseg(2:end-3);
    %x2  = xseg(4:end-1);
    %y2  = yseg(4:end-1);
    x = xseg(4:end-3);
    y = yseg(4:end-3);
    x1  = xseg(3:end-4);
    y1  = yseg(3:end-4);
    x2  = xseg(5:end-2);
    y2  = yseg(5:end-2);
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

coord_int = struct;
coord_ext = struct;
coord_int.x = x_int;
coord_int.y = y_int;
coord_ext.x = x_ext;
coord_ext.y = y_ext;
segment.coord_int = coord_int;
segment.coord_ext = coord_ext;
segment.ymax = max(segy);