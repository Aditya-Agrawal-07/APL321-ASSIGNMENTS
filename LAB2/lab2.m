% Question 1

Pe = 50;                      % Setting Peclet No. as 50
phi_0 = 0;                    % Boundary Condition at x = 0
phi_L = 1;                    % Boundary Condition at x = L

N = [11, 21, 41, 81];         % Set of total no. of grid points

uds = cell(1,size(N,2));      % array of arrays to store solutions of UDS for convective term
cds = cell(1,size(N,2));      % array of arrays to store solutions of CDS for convective term

for i = 1:size(N,2)
    
    n = N(1,i);               % No. of grid points
    h = 1/(n-1);              % step size

    % UDS solution assuming u > 0
    % Setting the parameters of UDS
    Ae = 1;                   % coeff. of phi(i+1)
    Ap = -(2 + Pe*h);         % coeff. of phi(i)
    Aw = 1 + Pe*h;            % coeff. of phi(i-1)

    % Defining the required sparse matrix A of size nxn
    A = zeros(n,n);
    A(1,1) = 1; A(n,n) = 1;
    for j = 2:n-1
        A(j,j) = Ap;
        A(j,j-1) = Aw;
        A(j,j+1) = Ae;
    end

    % Defining the load vector b of size nx1
    b = zeros(n,1);
    b(1,1) = phi_0; b(n,1) = phi_L;
    
    phi = linsolve(A, b);   % Using the inbuilt fnc linsolve to solve system of equations
    uds{1,i} = phi;         % Storing the solution in array of arrays


    % CDS solution
    % Setting the parameters of CDS
    Ae = 1 - (Pe*h)/2;      % coeff. of phi(i+1)                 
    Ap = -2;                % coeff. of phi(i)
    Aw = 1 + (Pe*h)/2;      % coeff. of phi(i-1)

    % Defining the required sparse matrix A of size nxn
    A = zeros(n,n);
    A(1,1) = 1; A(n,n) = 1;
    for j = 2:n-1
        A(j,j) = Ap;
        A(j,j-1) = Aw;
        A(j,j+1) = Ae;
    end

    % Defining the load vector b of size nx1
    b = zeros(n,1);
    b(1) = phi_0; b(n) = phi_L;
    
    phi = linsolve(A, b);   % Using the inbuilt fnc linsolve to solve system of equations
    cds{1,i} = phi;         % Storing the solution in array of arrays
end   

x = linspace(0,1,100);
phi_exact = (exp(x*Pe)-1)/(exp(Pe)-1);  % exact solution of the convection-diffusion eqn

% Plotting exact and computed solutions using UDS
figure;
plot(x,phi_exact,'r-','LineWidth',2)
for i=1:size(N,2)
    hold on
    interval = linspace(0,1,N(1,i));
    plot(interval, uds{1,i}, 'LineWidth', 1);
end

title('ɸ-ɸ_0/ɸ_L-ɸ_0 Vs x/L using UDS');
xlabel('x/L');
ylabel('ɸ-ɸ_0/ɸ_L-ɸ_0');
legend('Exact soln','N = 11','N = 21','N = 41', 'N = 81');
grid on;
xlim([0.000 1.000]);
ylim([-0.00 1.00]);
legend("Position", [0.15656,0.68826,0.20763,0.20984]);


% Plotting exact and computed solutions using CDS
figure;
plot(x,phi_exact,'r-','LineWidth',2)
for i=1:size(N,2)
    hold on
    interval = linspace(0,1,N(1,i));
    plot(interval, cds{1,i},'LineWidth',1);
end

title('ɸ-ɸ_0/ɸ_L-ɸ_0 Vs x/L using CDS');
xlabel('x/L');
ylabel('ɸ-ɸ_0/ɸ_L-ɸ_0');
legend('Exact soln','N = 11','N = 21','N = 41', 'N = 81');
legend('show');
grid on;
xlim([0.000 1.000]);
ylim([-0.50 1.00]);
legend("Position", [0.15497,0.70505,0.20634,0.19689]);

% Question 2

h_vals = zeros(size(N,2),1);        % array to store the h values
err_uds = zeros(size(N,2),1);       % array to store the discretization error of UDS
err_cds = zeros(size(N,2),1);       % array to store the discretization error of CDS

for i=1:size(N,2)
    n = N(1,i);                     % No. of grid points
    h = 1/(n-1);                    % step size
    h_vals(i) = h;
    x = linspace(0,1,n);
    x = transpose(x);
    phi_exact = (exp(x*Pe)-1)/(exp(Pe)-1);
    err_uds(i) = norm(uds{1,i}-phi_exact)/sqrt(n);
    err_cds(i) = norm(cds{1,i}-phi_exact)/sqrt(n);
end

c1 = 1.2;                           % coeff. of linear reference line
c2 = 12;                            % coeff. of quadratic reference curve

% Plotting the required graphs
figure;
loglog(h_vals, err_uds, 'LineWidth', 2);
hold on 
loglog(h_vals,err_cds,"LineWidth", 2);
hold on
plot(h_vals,c1*h_vals,'--', "LineWidth", 0.5);
hold on 
plot(h_vals,c2*h_vals.*h_vals, '--', "LineWidth",0.5);
xlabel('h');
ylabel('∈_h');
title('∈_h Vs h');
legend('UDS', 'CDS', 'y=c_1h', 'y=c_2h^2');
grid on;
legend("Position", [0.71467,0.15712,0.163,0.17487]);

% Question 3

% Computing the solutions for N = 161 and 321 as well
N = [11, 21, 41, 81, 161, 321];

uds = cell(1,size(N,2));      % array of arrays to store solutions of UDS for convective term
cds = cell(1,size(N,2));      % array of arrays to store solutions of CDS for convective term

for i = 1:size(N,2)
    
    n = N(1,i);               % No. of grid points
    h = 1/(n-1);              % step size

    % UDS solution assuming u > 0
    % Setting the parameters of UDS
    Ae = 1;                   % coeff. of phi(i+1)
    Ap = -(2 + Pe*h);         % coeff. of phi(i)
    Aw = 1 + Pe*h;            % coeff. of phi(i-1)

    % Defining the required sparse matrix A of size nxn
    A = zeros(n,n);
    A(1,1) = 1; A(n,n) = 1;
    for j = 2:n-1
        A(j,j) = Ap;
        A(j,j-1) = Aw;
        A(j,j+1) = Ae;
    end

    % Defining the load vector b of size nx1
    b = zeros(n,1);
    b(1,1) = phi_0; b(n,1) = phi_L;
    
    phi = linsolve(A, b);   % Using the inbuilt fnc linsolve to solve system of equations
    uds{1,i} = phi;         % Storing the solution in array of arrays


    % CDS solution
    % Setting the parameters of CDS
    Ae = 1 - (Pe*h)/2;      % coeff. of phi(i+1)                 
    Ap = -2;                % coeff. of phi(i)
    Aw = 1 + (Pe*h)/2;      % coeff. of phi(i-1)

    % Defining the required sparse matrix A of size nxn
    A = zeros(n,n);
    A(1,1) = 1; A(n,n) = 1;
    for j = 2:n-1
        A(j,j) = Ap;
        A(j,j-1) = Aw;
        A(j,j+1) = Ae;
    end

    % Defining the load vector b of size nx1
    b = zeros(n,1);
    b(1) = phi_0; b(n) = phi_L;
    
    phi = linsolve(A, b);   % Using the inbuilt fnc linsolve to solve system of equations
    cds{1,i} = phi;         % Storing the solution in array of arrays
end 

m_uds = zeros(4,1);         % array to store the order of accuracy corresponding to UDS
m_cds = zeros(4,1);         % array to store the order of accuracy corresponding to CDS
h_vals = zeros(4,1);        % array to store corresponding h values

for i = 3:6
    n1 = N(1,i);            % N corresponding to h
    n2 = N(1,i-1);          % N corresponding to 2h
    n4 = N(1,i-2);          % N corresponding to 4h

    h = 1/(n1-1);           % step size
    h_vals(i-2,1) = h;

    % computing the order of accuracy m_h for uds
    num = (norm(uds{1,i-1}(1:2:end) - uds{1,i-2}))/sqrt(n4);
    denom = (norm(uds{1,i}(1:2:end) - uds{1,i-1}))/sqrt(n2);
    m_uds(i-2,1) = log(num/denom)/log(2);

    % computing the order of accuracy m_h for cds
    num = (norm(cds{1,i-1}(1:2:end) - cds{1,i-2}))/sqrt(n4);
    denom = (norm(cds{1,i}(1:2:end) - cds{1,i-1}))/sqrt(n2);
    m_cds(i-2,1) = log(num/denom)/log(2);
end

% Plotting the m_h vs h curves for UDS and CDS
figure;
plot(h_vals,m_uds,'LineWidth',1);
hold on
plot(h_vals,m_cds,'LineWidth',1);
title('Order of Accuracy (m_h) Vs h');
legend('UDS','CDS');
xlabel('h');
ylabel('m_h');
grid on;