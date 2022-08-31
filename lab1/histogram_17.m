clc;
clear;
close all;

% h = csvread('histY00.csv');
% stairs(h(:,1),h(:,2), 'Color', 'k');
% title('Y00');
% grid on;
% figure;
% 
% h = csvread('histY01.csv');
% stairs(h(:,1),h(:,2), 'Color', 'k');
% title('Y01');
% grid on;
% figure;
% 
% h = csvread('histY10.csv');
% stairs(h(:,1),h(:,2), 'Color', 'k');
% title('Y10');
% grid on;
% figure;
% 
% h = csvread('histY11.csv');
% stairs(h(:,1),h(:,2), 'Color', 'k');
% title('Y11');
% grid on;

h = csvread('histD_Y00.csv');
stairs(h(:,1),h(:,2), 'Color', 'k');
title('D_Y_0_0');
grid on;
figure;

h = csvread('histD_Y01.csv');
stairs(h(:,1),h(:,2), 'Color', 'k');
title('D_Y_0_1');
grid on;
figure;

h = csvread('histD_Y10.csv');
stairs(h(:,1),h(:,2), 'Color', 'k');
title('D_Y_1_0');
grid on;
figure;

h = csvread('histD_Y11.csv');
stairs(h(:,1),h(:,2), 'Color', 'k');
title('D_Y_1_1');
grid on;