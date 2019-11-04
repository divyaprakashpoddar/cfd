% MATLAB Code for soving 1D Diffusion equation
% Author: Divyaprakash
% Course: Computational Fluid Dynamics

% Given data
gamma       = 0.5;          % Diffusion Coefficient
area_cv     = 1;            % Cross-sectional area
q           = 1000e3;

% Step 1: Create a grid
L           = 0.02;               % Length of domain (in meters)
N           = 25;                 % Number of grid points. Minimum value of N = 2

% Set the initial temperature distribution and Apply boundary conditions
T           = zeros(1,N+2);       % Temperature distribution
T(1)        = 100;             % Boundary Conditions
T(end)      = 200;

[myGrid, delX] = create_grid(L,N);
source      = [q*area_cv*delX 0];         % source = [Su Sp]; S = Su + Sp*phi_p
%source     = [0 0];         % source = [Su Sp]; S = Su + Sp*phi_p


% Step 2: Discretize the equations at each control volue
% around a node and form the coefficient matrix
A           = zeros(N);           % Coefficient Matrix
b           = zeros(N,1);         % RHS Vector

% Fill up matrix A; 
% -Aw*Tw + Ap*Tp -Ae*Te = Su

for P = 2:N-1
    [aw, ap, ae, su] = get_coefficients(N, delX, P, gamma, area_cv, source, T);
    W = P - 1;
    E = P + 1;
    A(P,W:E) = [-aw, ap, -ae];
    b(P) = su;
end

% Boundaries
[aw, ap, ae, su] = get_coefficients(N, delX, 1, gamma, area_cv, source, T);
A(1,1:2) = [ap, -ae];
b(1) = su;
[aw, ap, ae, su] = get_coefficients(N, delX, N, gamma, area_cv, source, T);
A(N,N-1:N) = [-aw, ap];
b(N) = su;


% Step 3: Solution of equations. Solve for the unknowns; Ax = b
%Temp = inv(A)*b;
Temp = A\b;

T(2:end-1) = Temp;

% Exact Solution
X = 0:0.001:L;
T_exact = ((T(end) - T(1))/L + q/(2*gamma)*(L - X)).*X + T(1);

hold all
grid on
plot(X, T_exact, 'g');
scatter(myGrid, T);
legend('Exact', 'Numerical','location','northwest');
xlabel('Distance x(m)');
ylabel('Temperature ({}^0{C})')

