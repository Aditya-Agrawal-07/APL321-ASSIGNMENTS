% Discretization of convection-diffusion equation using CDS scheme

phi_0 = 0;
phi_L = 1;
Pe = 50;                                % Peclet No. = 50
N = 41;                                 % Grid Points = 41
h = 1/(N-1);                            % Grid Width = 1/(N-1)

% Setting the parameters of CDS
Ae = 1 - (Pe*h)/2;      % coeff. of phi(i+1)                 
Ap = -2;                % coeff. of phi(i)
Aw = 1 + (Pe*h)/2;      % coeff. of phi(i-1)

% Defining the required sparse matrix A of size nxn
A = zeros(N,N);
A(1,1) = 1; A(N,N) = 1;
for j = 2:N-1
    A(j,j) = Ap;
    A(j,j-1) = Aw;
    A(j,j+1) = Ae;
end

% Defining the load vector b of size nx1
b = zeros(N,1);
b(1) = phi_0; b(N) = phi_L;

% Using the inbuilt fnc linsolve to solve system of equations
phi_h = linsolve(A, b);           

x = transpose(linspace(0,1,41));
% exact solution of the convection-diffusion eqn
phi_exact = (exp(x*Pe)-1)/(exp(Pe)-1);

% Question 2

Pe = 50;                                % Peclet No. = 50
N = 41;                                 % Grid Points = 41
h = 1/(N-1);                            % Grid Width = 1/(N-1)
err = zeros(N,1);                       % vector to store truncation error
x = transpose(linspace(0,1,N));         % x_i values
phi3 = (Pe^3*exp(x*Pe))/(exp(Pe)-1);    % third order derivative values of phi_exact
phi4 = (Pe^4*exp(x*Pe))/(exp(Pe)-1);    % fourth order derivative values of phi_exact

for i = 2:N-1
    err(i) = (Pe*h*h/6)*phi3(i) - (h*h/12)*phi4(i); % substituting the value of err(i)
end

% Plot of truncation error vs x_i
figure;
plot(x, err, 'LineWidth',1);
title('Truncation error(∈_τ) Vs x/L')
xlabel('x/L');
ylabel('∈_τ');
grid on;

% Question 3

err_h = phi_exact - phi_h;      % discretization error

L = zeros(N,N);
L(1,1) = 1; L(N,N) = 1;
Aw = 1/h^2 + Pe/(2*h);
Ae = 1/h^2 - Pe/(2*h);
Ap = -2/h^2;

for j = 2:N-1
    L(j,j) = Ap;
    L(j,j-1) = Aw;
    L(j,j+1) = Ae;
end

err_th = -inv(L)*err;           % approx. discretization error

% Plotting exact and approx. discretization errors vs h

figure;
plot(x,err_h,'LineWidth',1);
hold on
plot(x, err_th,'LineWidth',1);
title('Exact and Approx. Discretization Errors Vs x/L')
xlabel('x/L');
ylabel('Error (∈)');
legend('Exact error (∈_h)','Approx. error (∈_{th})');
xlim([0.00 1.00])
ylim([-0.0100 0.0600])
legend("Position", [0.18496,0.7564,0.3824,0.12664]);
grid on;

% Question 4

% Solving for phi_exact
N = [41; 81; 161; 321];
phi_exact = cell(size(N,1),1);

for i=1:size(N,1)
    n = N(i);
    x = transpose(linspace(0,1,n));
    phi = (exp(x*Pe)-1)/(exp(Pe)-1);
    phi_exact{i,1} = phi;
end

% Solving for phi_h
phi_h = cell(size(N,1),1);

for i=1:size(N,1)
    n = N(i);                                 % Grid Points = 41
    h = 1/(n-1);                            % Grid Width = 1/(N-1)
    
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
    
    % Using the inbuilt fnc linsolve to solve system of equations
    soln = linsolve(A, b);
    phi_h{i,1} = soln;
end

% Solving for err_h

err_h = cell(size(N,1),1);

for i = 1:size(N,1)
    err_h{i,1} = phi_exact{i,1} - phi_h{i,1};
end

% Solving for the std. dev. norm of err_h
err_h_norm = zeros(size(N,1),1);

for i = 1:size(N,1)
    val = norm(err_h{i,1})/sqrt(N(i));
    err_h_norm(i) = val;
end

% Solving for truncation_error

trunc_error = cell(size(N,1),1);

for i = 1:size(N,1)
    n = N(i);                               % Grid Points = 41
    h = 1/(n-1);                            % Grid Width = 1/(N-1)
    x = transpose(linspace(0,1,n));
    phi3 = (Pe^3*exp(x*Pe))/(exp(Pe)-1);    % third order derivative values of phi_exact
    phi4 = (Pe^4*exp(x*Pe))/(exp(Pe)-1);    % fourth order derivative values of phi_exact
    vec = zeros(n,1);
    for j = 2:n-1
        vec(j) = (Pe*h*h/6)*phi3(j) - (h*h/12)*phi4(j); % substituting the value of err(i)
    end
    trunc_error{i,1} = vec;
end

% Solving for err_th

err_th = cell(size(N,1),1);
h_vals = zeros(size(N,1),1);

for i = 1:size(N,1)
    n = N(i);                               % Grid Points = 41
    h = 1/(n-1);                            % Grid Width = 1/(N-1)
    h_vals(i) = h;
    L = zeros(n,n);
    L(1,1) = 1; L(n,n) = 1;
    Aw = 1/h^2 + Pe/(2*h);
    Ae = 1/h^2 - Pe/(2*h);
    Ap = -2/h^2;
    
    for j = 2:n-1
        L(j,j) = Ap;
        L(j,j-1) = Aw;
        L(j,j+1) = Ae;
    end
    
    err_th{i,1} = -inv(L)*trunc_error{i,1};           % approx. discretization error
end

% Solving for the std. dev. norm of err_th
err_th_norm = zeros(size(N,1),1);

for i = 1:size(N,1)
    val = norm(err_th{i,1})/sqrt(N(i));
    err_th_norm(i) = val;
end

% Plotting err_h_norm and err_th_norm vs x_i
c = 12;

figure;
loglog(h_vals, err_h_norm, 'Linewidth', 2);
hold on
loglog(h_vals, err_th_norm, 'Linewidth', 2);
hold on
plot(h_vals,c*h_vals.*h_vals, '--', "LineWidth",0.5);
xlabel('h');
ylabel('||∈_h||');
title('Std. Dev. Norm of Discretization Errors Vs h');
legend('Exact Error', 'Approx. Error', 'y=ch^2');
grid on;
legend("Position", [0.57785,0.16243,0.30196,0.17396])
