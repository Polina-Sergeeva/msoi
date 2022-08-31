clc;
clear;
close all;

h = csvread('histR.csv');
stairs(h(:,1),h(:,2), 'Color', 'k');
title('R');
grid on;
figure;

h = csvread('histG.csv');
stairs(h(:,1),h(:,2), 'Color', 'k');
title('G');
grid on;
figure;

h = csvread('histB.csv');
stairs(h(:,1),h(:,2), 'Color', 'k');
title('B');
grid on;
figure;

h = csvread('histY.csv');
stairs(h(:,1),h(:,2), 'Color', 'k');
title('Y');
grid on;
figure;

h = csvread('histCb.csv');
stairs(h(:,1),h(:,2), 'Color', 'k');
title('Cb');
grid on;
figure;

h = csvread('histCr.csv');
stairs(h(:,1),h(:,2), 'Color', 'k');
title('Cr');
grid on;