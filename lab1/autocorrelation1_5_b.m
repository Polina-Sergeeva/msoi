for m = 1:5
    comp = csvread(['autocorrelationY_', num2str(m), '.csv']);
    plot3(comp(:,1), comp(:,2), comp(:,3));
    hold on;
end
grid on;
xlabel('x');
ylabel('y');
zlabel('r(x,y)'); 
figure;

for m = 1:5
    comp = csvread(['autocorrelationCb_', num2str(m), '.csv']);
    plot3(comp(:,1), comp(:,2), comp(:,3));
    hold on;
end
grid on;
xlabel('x');
ylabel('y');
zlabel('r(x,y)'); 
figure;

for m = 1:5
    comp = csvread(['autocorrelationCr_', num2str(m), '.csv']);
    plot3(comp(:,1), comp(:,2), comp(:,3));
    hold on;
end
grid on;
xlabel('x');
ylabel('y');
zlabel('r(x,y)') ;