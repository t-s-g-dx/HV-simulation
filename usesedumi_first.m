% Define variables
% alpha1 = sdpvar(1, 1);
alpha1=1;
alpha2 = sdpvar(1, 1);
% alpha2=1;
alpha3 = sdpvar(1, 1);
% alpha3=1;
k1 = 3;
k2 = 5;
% k1=sdpvar(1,1);
% k2=sdpvar(1,1);
% Define small positive constant
epsilon = 1e-20;


% Define LMIs with non-strict inequalities
LMI1 = 3*alpha1*k1-2*alpha2*k2 >= epsilon;

LMI2 = alpha2-3*alpha3*k2 >= epsilon;
LMI3 = alpha3-alpha2 >= epsilon;
LMI4 = alpha1-2*alpha2 >= epsilon;
LMI5 = 4*alpha2+12*alpha3*k2-3*alpha1-2*alpha2*k2 >= epsilon;

LMI6 = 3*alpha1*k1+6*alpha3*k2+2*alpha2-6*alpha1-6*alpha2*k2>=epsilon;

% Combine LMIs
LMIs = [LMI1, LMI2, LMI3, LMI4, LMI5,LMI6];
% Set optimization options
options = sdpsettings('solver', 'sedumi', 'verbose', 0);

% Solve the LMI problem
diagnostics = optimize(LMIs, [], options);

if diagnostics.problem == 0
    disp('Optimization successful!');
    disp('Optimal solution:');
    disp(['alpha1 = ', num2str(value(alpha1))]);
    disp(['alpha2 = ', num2str(value(alpha2))]);
    disp(['alpha3 = ', num2str(value(alpha3))]);
else
    disp('Optimization failed or no optimal solution found!');
end