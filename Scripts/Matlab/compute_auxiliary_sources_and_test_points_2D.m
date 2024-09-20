function segment = compute_auxiliary_sources_and_test_points_2D(segx,segy,alpha)

d = segx(2)-segx(1);
r = min((max(segx) - min(segx)),(max(segy) - min(segy)))/2;

ys = 0:d:segy(end);
x = [segx 0*ys+max(segx)];
y = [segy flip(ys)];
x = [x flip(segx)];
y = [y 0*segy];
ys = 0:d:segy(1);
x = [x 0*ys+min(segx)];
y = [y ys];

x_left  = [x(end) x(1:end-1)];
y_left  = [y(end) y(1:end-1)];
x_right = [x(2:end) x(1)];
y_right = [y(2:end) y(1)];
x_diff = x_left - x_right;
y_diff = y_left - y_right;
diff = [x_diff; y_diff];
diff = diff./vecnorm(diff)*r*(1-alpha);

x_int = x - diff(2,:);
y_int = y + diff(1,:);
x_ext = x + diff(2,:);
y_ext = y - diff(1,:);

% Save data in struct
segment = struct;
segment.x = x;
segment.y = y;
segment.x_int = x_int;
segment.y_int = y_int;
segment.x_ext = x_ext;
segment.y_ext = y_ext;