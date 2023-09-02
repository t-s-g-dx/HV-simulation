% Define variables
gamma12 = sdpvar(1, 1);
gamma23 = sdpvar(1, 1);
gamma1 = sdpvar(1, 1);
gamma2 = sdpvar(1, 1);
gamma3 = sdpvar(1, 1);
% gamma13 =sdpvar(1, 1);
gamma13=0;
k1 = 2;
k2 = 1;
k3 = 1;
% k1=sdpvar(1,1);
% k2=sdpvar(1,1);
% k3=sdpvar(1,1);
% Define small positive constant
epsilon = 1e-7;

% Build the matrix
Gamma = [gamma1, -1/2 * gamma12, 0; -1/2 * gamma12, gamma2, -1/2 * gamma23; 0, -1/2 * gamma23, gamma3];

% Define LMIs with non-strict inequalities
LMI1 = Gamma >= epsilon * eye(3);
LMI2 = gamma12 >= epsilon;
LMI3 = 4/3 * gamma1 * k1 - 2/3 * gamma12 * k2 - 4/9 * gamma1 - 2/9 * gamma12 * k1 + 2/3 * gamma13 >= epsilon;
LMI4 = 2 * gamma2 * k2 - 1/3 * gamma12 * k2 - 8/9 * gamma1 - 4/9 * gamma12 * k1 - 1/3 * gamma23 * k2 - 4/3 * gamma2 - 4/3 * gamma23 * k3 >= epsilon;
LMI5 = 4 * gamma3 * k3 + gamma13 - 2/3 * gamma23 * k2 - 2/3 * gamma2 - 2/3 * gamma23 * k3 >= epsilon;

% Combine LMIs
LMIs = [LMI1, LMI2, LMI3, LMI4, LMI5];

% Set optimization options
options = sdpsettings('solver', 'sedumi', 'verbose', 0);

% Solve the LMI problem
diagnostics = optimize(LMIs, [], options);

if diagnostics.problem == 0
    disp('Optimization successful!');
    disp('Optimal solution:');
    disp(['gamma12 = ', num2str(value(gamma12))]);
    disp(['gamma23 = ', num2str(value(gamma23))]);
    disp(['gamma1 = ', num2str(value(gamma1))]);
    disp(['gamma2 = ', num2str(value(gamma2))]);
    disp(['gamma3 = ', num2str(value(gamma3))]);
    disp(['gamma13 = ', num2str(value(gamma13))]);
else
    disp('Optimization failed or no optimal solution found!');
end