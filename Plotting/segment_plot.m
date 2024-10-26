clear; close all; clc;


num_segments = 10;

figure('Renderer', 'painters', 'Position', [400 400 1500 500]);

for i = 1:1
    
    data = readmatrix(sprintf('../Data/segments/test_top_segment_%d.txt',i));
    x = data(:,1);
    y = data(:,2);

    plot(x,y,'.')
    hold on

    data = readmatrix(sprintf('../Data/segments/test_n_segment_%d.txt',i));
    nx = data(:,1);
    ny = data(:,2);

    plot(x+nx*10^(-9),y+ny*10^(-9),'ko','MarkerSize',10,'LineWidth',1)
    hold on

    data = readmatrix(sprintf('../Data/segments/test_right_segment_%d.txt',i));
    x = data(:,1);
    y = data(:,2);

    plot(x,y,'.')
    hold on

    data = readmatrix(sprintf('../Data/segments/test_bottom_segment_%d.txt',i));
    x = data(:,1);
    y = data(:,2);

    plot(x,y,'.')
    hold on

    data = readmatrix(sprintf('../Data/segments/test_left_segment_%d.txt',i));
    x = data(:,1);
    y = data(:,2);

    plot(x,y,'o')
    hold on

    data = readmatrix(sprintf('../Data/segments/test_int_segment_%d.txt',i));
    x = data(:,1);
    y = data(:,2);

    plot(x,y,'*')
    hold on

    

    
    axis equal

    xlim([0,0.25]*10^(-6))
    ylim([-1,4]*10^(-8))

end