clc;
clear;
close all;

h = csvread('histD_R.csv');
stairs(h(:,1),h(:,2), 'Color', 'k');
title('D_R');
grid on;
figure;

h = csvread('histD_G.csv');
stairs(h(:,1),h(:,2), 'Color', 'k');
title('D_G');
grid on;
figure;

h = csvread('histD_B.csv');
stairs(h(:,1),h(:,2), 'Color', 'k');
title('D_B');
grid on;
figure;

h = csvread('histD_Y.csv');
stairs(h(:,1),h(:,2), 'Color', 'k');
title('D_Y');
grid on;
figure;

h = csvread('histD_Cb.csv');
stairs(h(:,1),h(:,2), 'Color', 'k');
title('D_C_b');
grid on;
figure;

h = csvread('histD_Cr.csv');
stairs(h(:,1),h(:,2), 'Color', 'k');
title('D_C_r');
grid on;