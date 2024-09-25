
clear; close all; clc;

filename = "demoleus2x2"; % Retinin2x2 demoleus2x2
views = {"far","close"};
view = "close"; % close far

total = 1000;

for hn = 1:2
view = views{hn};
for num_segments = [1 5 10 20]
load(sprintf("../../Data/%s_2D.mat",filename))
m = length(X);
%num_segments = 20; % 1 5 10 20
n = round(m/num_segments);
n = round(total/num_segments);
m = num_segments * n;

x = linspace(min(X),max(X),m+1);
d = x(2)-x(1);

r = (max(x) - min(x))/(2*num_segments);
alpha = 2*d;

y = interp1(X,Y,x);

figure('Renderer', 'painters', 'Position', [400 400 800 250]);
plot(x,y,'LineWidth',1.5)
xlim([min(x),max(x)])
grid on
hold on

segments = cell(num_segments,1);

figure('Renderer', 'painters', 'Position', [400 400 1000 800]);
for k = 1:num_segments
    
    segx = x((k-1)*n+1:k*n+1);
    segy = y((k-1)*n+1:k*n+1);
    
    segments{k} = setup_nanostructures_2D(segx,segy,alpha);
    for j =1:length(segments{k}.x)
        plot(segments{k}.x{j},segments{k}.y{j},'r.','LineWidth',1.5,'Markersize',10)
        hold on
    end
    
    plot(segments{k}.x_int,segments{k}.y_int,'b.','LineWidth',1.5,'Markersize',10)
    
    plot(segments{k}.x_ext,segments{k}.y_ext,'go','LineWidth',1.5,'Markersize',5)
    plot(segments{k}.n_x,segments{k}.n_y,'k.','LineWidth',1.5,'Markersize',10)
end
xlim([-0.5*10^(-7),20.5*10^(-7)])
ylim([-0.5*10^(-7),5*10^(-7)])
grid on
axis equal

tic;
for k = 1:length(segments)
    segments{k} = forward_nanostructures_2D(segments{k});
end
toc;
n = 300;


X = linspace(0,21*10^(-7),n);
Y = linspace(-0.5*10^(-7),20.5*10^(-7),n);
if strcmp("far",view)
    Y = Y + 3*10^(-2);
end
[Xmesh,Ymesh] = meshgrid(X,Y);


[Etot, Htot, Einc, Hinc, Eref, Href, Escat, Hscat] = compute_2D_field_nanostructures_2D(Xmesh,Ymesh,segments);

filenamesave = sprintf('../../Results/nanostructures_2D/%s_%s_num_segments_%d.mat',filename,view,num_segments);
save(filenamesave,'Etot','Htot','Einc','Hinc','Eref','Href','Escat','Hscat','X','Y','segments');
end
end
%%
figure;
imagesc(x,y,abs(Etot(:,:,3)).^2)
set(gca,'YDir','normal')
colorbar
hold on
for j = 1:length(segments{1}.x)
    %plot(segments{1}.x{j},segments{1}.y{j},'r.','LineWidth',1.5,'Markersize',10)
end
axis equal
axis tight