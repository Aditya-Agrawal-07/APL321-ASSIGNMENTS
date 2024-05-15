close all;
clear all;
clc;

% Question 3

% Defining the constants of the equation

Lx = 20; Ly = 1;        % dimensions of rectangle
D = 1e-4;               % species diffusivity
u0 = 0.01;              % max channel velocity
nx = 401; ny = 21;      % no. of grid points
x = linspace(0,Lx,nx);  % x coordinates of nodes
y = linspace(0,Ly,ny);  % y coordinates of nodes
xc = x(101);            % location of initial concentration
w = [1.05; 1.05];       % omega for SOR method
U = f(y, u0, Ly);       % velocity field


% Scheme 1
c1 = zeros(ny,nx);      % solution matrix
dx = x(2) - x(1); dy = y(2) - y(1); % grid width
c1(:,100:102) = 1/(2*dx);           % initial conc.

dt = 2;                 % time step

% Coefficients of discretized eqn
Ae = -D/(dx^2); An = -D/(dy^2); As = -D/(dy^2); 
Aw = -(U/dx + D/(dx^2)); 
Ap = 1/dt + U/dx + 2*D/(dx^2) + 2*D/(dy^2);

for t=1:500/dt
    
    c_prev = c1;
    err = 100;
    w1 = w(1);
    iter = 0;

    while (abs(err) > 1e-10)

        c0 = c1;

        for i=2:nx-1
            for j=2:ny-1
                Q = c_prev(j,i)/dt;
                c1(j,i) = (Q - (Ae*c1(j,i+1) + Aw(j)*c1(j,i-1) + An*c1(j+1,i) + As*c1(j-1,i)))/Ap(j);
                d = c1(j,i) - c0(j,i);
                c1(j,i) = c0(j,i) + w1*d;
       
            end
        end

        c1(1,:) = c1(2,:);
        c1(:,nx) = c1(:,nx-1);
        c1(ny,:) = c1(ny-1,:);

        err = norm(c1 - c0)/sqrt(nx*ny);
        iter = iter + 1;
    end
end

% plotting the isocontours of solution
figure;
contourf(x,y,c1);
title("Isocontours for c(x,y,500) using Scheme 1");
xlabel("x");
ylabel("y");
colorbar;



% Scheme 2
c2 = zeros(ny,nx);      % solution matrix
dx = x(2) - x(1); dy = y(2) - y(1); % grid width
c2(:,100:102) = 1/(2*dx);           % initial conc.

dt = 2;                 % time step

% Coefficients of discretized eqn
Ae = -D/(2*dx^2); An = -D/(2*dy^2); As = -D/(2*dy^2);
Aw = -(U/dx + D/(2*dx^2)); 
Ap = 1/dt + U/dx + D/(dx^2) + D/(dy^2);


for t=1:500/dt
    
    c_prev = c2;
    err = 100;
    w1 = w(2);
    iter = 0;

    while (abs(err) > 1e-10)

        c0 = c2;

        for i=2:nx-1
            for j=2:ny-1
                Q = (1/dt - D/(dx^2) - D/(dy^2))*c_prev(j,i) + (D/2)*((c_prev(j,i+1)+c_prev(j,i-1))/(dx^2) + (c_prev(j+1,i)+c_prev(j-1,i))/(dy^2));
                c2(j,i) = (Q - (Ae*c2(j,i+1) + Aw(j)*c2(j,i-1) + An*c2(j+1,i) + As*c2(j-1,i)))/Ap(j);
                d = c2(j,i) - c0(j,i);
                c2(j,i) = c0(j,i) + w1*d;
       
            end
        end

        c2(1,:) = c2(2,:);
        c2(:,nx) = c2(:,nx-1);
        c2(ny,:) = c2(ny-1,:);

        err = norm(c2 - c0)/sqrt(nx*ny);
        iter = iter + 1;
    end
end

% plotting the isocontours of solution
figure;
contourf(x,y,c2);
title("Isocontours for c(x,y,500) using Scheme 2");
xlabel("x");
ylabel("y");
colorbar;




% Question 4

delta_t = [0.5; 1; 1.5; 2];
scheme1 = cell(4,1); scheme2 = cell(4,1);

for k=1:4

    % Scheme 1
    c1 = zeros(ny,nx);
    dx = x(2) - x(1); dy = y(2) - y(1);
    c1(:,100:102) = 1/(2*dx);
    
    dt = delta_t(k);
    
    Ae = -D/(dx^2); An = -D/(dy^2); As = -D/(dy^2);
    Aw = -(U/dx + D/(dx^2)); 
    Ap = 1/dt + U/dx + 2*D/(dx^2) + 2*D/(dy^2);
    
    
    for t=1:500/dt
        
        c_prev = c1;
        err = 100;
        w1 = w(1);
        iter = 0;
    
        while (abs(err) > 1e-10)
    
            c0 = c1;
    
            for i=2:nx-1
                for j=2:ny-1
                    Q = c_prev(j,i)/dt;
                    c1(j,i) = (Q - (Ae*c1(j,i+1) + Aw(j)*c1(j,i-1) + An*c1(j+1,i) + As*c1(j-1,i)))/Ap(j);
                    d = c1(j,i) - c0(j,i);
                    c1(j,i) = c0(j,i) + w1*d;
           
                end
            end
    
            c1(1,:) = c1(2,:);
            c1(:,nx) = c1(:,nx-1);
            c1(ny,:) = c1(ny-1,:);
    
            err = norm(c1 - c0)/sqrt(nx*ny);
            iter = iter + 1;
        end
    end

    scheme1{k,1} = c1(11,:);
    
    
    % Scheme 2
    c2 = zeros(ny,nx);
    dx = x(2) - x(1); dy = y(2) - y(1);
    c2(:,100:102) = 1/(2*dx);
    
    dt = delta_t(k);
    
    Ae = -D/(2*dx^2); An = -D/(2*dy^2); As = -D/(2*dy^2);
    Aw = -(U/dx + D/(2*dx^2)); 
    Ap = 1/dt + U/dx + D/(dx^2) + D/(dy^2);
    
    
    for t=1:500/dt
        
        c_prev = c2;
        err = 100;
        w1 = w(2);
        iter = 0;
    
        while (abs(err) > 1e-10)
    
            c0 = c2;
    
            for i=2:nx-1
                for j=2:ny-1
                    Q = (1/dt - D/(dx^2) - D/(dy^2))*c_prev(j,i) + (D/2)*((c_prev(j,i+1)+c_prev(j,i-1))/(dx^2) + (c_prev(j+1,i)+c_prev(j-1,i))/(dy^2));
                    c2(j,i) = (Q - (Ae*c2(j,i+1) + Aw(j)*c2(j,i-1) + An*c2(j+1,i) + As*c2(j-1,i)))/Ap(j);
                    d = c2(j,i) - c0(j,i);
                    c2(j,i) = c0(j,i) + w1*d;
           
                end
            end
    
            c2(1,:) = c2(2,:);
            c2(:,nx) = c2(:,nx-1);
            c2(ny,:) = c2(ny-1,:);
    
            err = norm(c2 - c0)/sqrt(nx*ny);
            iter = iter + 1;
        end
    end
        
    scheme2{k,1} = c2(11,:);

end

figure;
for k = 1:4
    plot(x,scheme1{k,1});
    if (k~=4)
        hold on;
    end
end
title("c(x,Ly/2,500) vs x for Scheme 1");
xlabel("x");
ylabel("c(x,Ly/2,500)");
legend("△t=0.5","△t=1","△t=1.5","△t=2");


figure;
for k = 1:4
    plot(x,scheme2{k,1});
    if (k~=4)
        hold on;
    end
end
title("c(x,Ly/2,500) vs x for Scheme 2");
xlabel("x");
ylabel("c(x,Ly/2,500)");
legend("△t=0.5","△t=1","△t=1.5","△t=2");

function U = f(y, u0, Ly)
    
    U = 4*u0*(y/Ly).*(1-y/Ly);
end
