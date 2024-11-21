clear; clc; close all; 

% Sample data
x = linspace(0, 2*pi, 100);
y1 = sin(x);
y2 = cos(x);
y3 = tan(x);
y4 = exp(-0.1*x).*sin(3*x);

f = figure(1); 
tiledlayout(f,1,4)
for i = 1:4 
    nexttile([i 2])
    subplot(1,2,1)
    plot(x,y1)
    subplot(1,2,2)
    plot(x,y2)
end