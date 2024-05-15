clear all;
close all;
clc;

% Question 5

load('solution.mat');       % loading the solution matrices obtained from Q4

y = [10; 40; 160];          % center of CV just below point y = 0.5 
x = cell(3,1); phi = cell(3,1); PHI = cell(3,1);
x{1,1} = x1; x{2,1} = x2; x{3,1} = x3;
phi{1,1} = phi1; phi{2,1} = phi2; phi{3,1} = phi3;
PHI{1,1} = PHI1; PHI{2,1} = PHI2; PHI{3,1} = PHI3;

% Plotting ɸ(x,y=0.5) vs x for UDS solution
figure;

for k=1:3

    plot(x{k,1}, (phi{k,1}(y(k),:)+phi{k,1}(y(k)+1,:))/2 , LineWidth=0.8);
    if k~= 3
        hold on;
    end
end

xlabel('x');
ylabel('ɸ');
title('ɸ Vs x for y = 0.5')
legend('N = 20', 'N = 80', 'N = 320');

% Plotting ɸ(x,y=0.5) vs x for CDS solution
figure;

for k=1:3

    plot(x{k,1}, (PHI{k,1}(y(k),:)+PHI{k,1}(y(k)+1,:))/2 , LineWidth=0.8);
    if k~= 3
        hold on;
    end
end

xlabel('x');
ylabel('ɸ');
title('ɸ Vs x for y = 0.5')
legend('N = 20', 'N = 80', 'N = 320');


