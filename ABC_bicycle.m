clc
clear all
%% Parameters
global par
par.L = 2.5; % Length of the bicycle
par.xi_0 = 20; % Horizontal coordinate of the obstacle
par.eta_0 = -0.1; % Vertical coordinate of the obstacle
par.R_0 = 4; % Radius of the obstacle
par.alpha_2 = 5; % Class-K for ABC
par.alpha = 1; % Class-K for virtual controller kappa

par.sigma = 0.001; % Smoothing parameter
par.mu = 1; 
par.t_max = 15; % Maximal time of the simulation
par.vd = 10; % Desired velocity
par.vd0 = 4; % Desired velocity in the desired virtual controller kappa
par.Keta = 0.4; % Gain in desired controller
par.Ktheta = 1.75; % Gain in desired controller
par.Kv = .3; % Gain in desired controller
%% Simulation
IC = [0, 0, 0, 2]; % Initial condition
opt = odeset("RelTol",1e-9,"AbsTol",1e-8);
[t,x] = ode45(@Bicycle_dynamics, [0 par.t_max], IC, opt);
%% Evaluation of functions
Safe_inputs = zeros(2, length(t));
kappa = zeros(2, length(t));
h_ABC = zeros(1, length(t));
psi_constr = zeros(1, length(t));
k_desired = zeros(2, length(t));
for ii=1:length(t)
    Safe_inputs(:,ii) = u_ABC(x(ii,(1:4))');
    kappa(:,ii) = k_0(x(ii,(1:4))');    
    h_ABC(:,ii) = ABC(x(ii,(1:4))');
    psi_constr(:,ii) = psi(x(ii,(1:4))');
    k_desired(:,ii) = k_d(x(ii,(1:4))'); 
end
%% First figure: Obstacle avoidance
figure(1)
axis equal
box on
grid off
hold on
f1 = @(x,y) (x-par.xi_0).^2 + (y-par.eta_0).^2 - par.R_0^2;
fp1=fimplicit(f1,[0.1,25,-4.5,4.5]);
hold on
fp2=fcontour(f1,[0.1,25,-4.5,4.5]);
fp2.LevelList = [0 0];
fp2.Fill = 'on';
colormap([[248 150 150]/256; 1 1 1])
hold on
fimplicit(@(x,y) (x-par.xi_0).^2 + (y-par.eta_0).^2 - par.R_0^2,'LineWidth',2,'Color','red')
hold on
plot(x(:,1), x(:,2),'LineWidth',2, 'Color', 'g')
ylim([-5 5])
xlim([0 35])
xlabel('$\xi\,$(m)','interpreter','latex')
ylabel('$\eta\,$(m)','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize', 10) 
plot(par.xi_0, par.eta_0, '.', 'LineWidth',8, 'Color', 'k')
yline(0)
%% Second figure: CBF
figure(2)
box on
grid off
f9 = @(x,y) -y;
fp9=fimplicit(f9,[0,15,-10,400]);
hold on
fp10=fcontour(f9,[0,15,-10,400]);
fp10.LevelList = [0 0];
fp10.Fill = 'on';
colormap([[185 217 190]/256; [185 217 190]/256; 1 1 1; [248 150 150]/256])
hold on
plot(t,h_ABC,'LineWidth',2,'Color','b')
hold on
plot(t,psi_constr,'LineWidth',2,'Color','r', 'LineStyle', ':')
hold on
xlabel('time, $t$ (s)','interpreter','latex')
ylabel('CBF, $h$ ($\mathrm{m^2}$)','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize', 10)
hold on
yline(0, 'LineWidth',1.5, 'Color', [0,170,0]/256)
ylim([-10,h_ABC(1)*1.01])
xlim([0 15])
legend('','','ABC', '$\psi$', 'interpreter', 'latex')
%% Third figure: State variables and inputs
figure(3)
subplot 221
box on
hold on
grid off
plot(t,Safe_inputs(1,:),'LineWidth',2,'Color','b')
hold on
plot(t,k_desired(1,:),'LineWidth',2,'Color','r','LineStyle','--')
hold on
ylabel('input ,$u_1\,$ (deg)','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize', 10) 
xlim([0 15])
legend('$u_1$', '$k_{\mathrm{d},1}$', 'interpreter', 'latex')

subplot 223
box on
hold on
grid off
plot(t,Safe_inputs(2,:),'LineWidth',2,'Color','b')
hold on
grid off
plot(t,k_desired(2,:),'LineWidth',2,'Color','r','LineStyle','--')
hold on
xlabel('time, $t$ (s)','interpreter','latex')
ylabel('input, $u_2\,$(m/$\mathrm{s^2})$','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize', 10) 
xlim([0 15])
legend('$u_2$', '$k_{\mathrm{d},2}$', 'interpreter', 'latex')

subplot 222
hold on
box on
grid off
plot(t,x(:,3),'LineWidth',2,'Color','b')
ylabel('yaw angle, $\vartheta$ (rad)','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize', 10) 
hold on
plot(t',zeros(length(t),1),'LineWidth',2,'Color','r','LineStyle','--')
legend("$\vartheta$","$\vartheta_{\mathrm{G}}$",'Interpreter','latex')
xlim([0 15])

subplot 224
hold on
box on
grid off
plot(t,x(:,4),'LineWidth',2,'Color','b')
xlabel('time, $t$ (s)','interpreter','latex')
ylabel('velocity, $v$ (m/s)','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize', 10) 
hold on
plot(t',zeros(length(t),1)+par.vd,'LineWidth',2,'Color','r','LineStyle','--')
legend("$v$","$v_{\mathrm{G}}$",'Interpreter','latex')
xlim([0 15])
%% Functions
function out = Bicycle_dynamics(t,x)
disp(t)
out = f(x) + g(x)*u_ABC(x);
end

function [Top,gradTop] = top_layer(x)
f1 = x(4)*cos(x(3));
f2 = x(4)*sin(x(3));
Top = [f1; f2];
gradTop = [0 0 -x(4)*sin(x(3)) cos(x(3));
           0 0  x(4)*cos(x(3)) sin(x(3))];
end

function out = f(x)
[Top,~] = top_layer(x);
out = [Top; 0; 0];
end

function out = g(x)
global par
out = [0 0;
       0 0;
       x(4)/par.L 0; 
       0 1];
end

function [Psi,gradPsi_top,gradPsi,gradgradPsi_top] = psi(x)
global par
Psi = (x(1)-par.xi_0)^2 + (x(2)-par.eta_0)^2 - par.R_0^2; % Constraint function
gradPsi_top = 2*[(x(1)-par.xi_0) (x(2)-par.eta_0)]; % Gradient with respect to variables in the top dynamics
gradPsi = [gradPsi_top 0 0]; % Gradient with respect to all variables
gradgradPsi_top = 2*[eye(2) zeros(2,2)]; % Hessian matrix with respect to variables in the top dynamics
end

function [K_d_0, gradK_d_0] = k_d_0(x) % Desired controller for the virtual controller kappa
global par
K_d_0 = [par.vd0;0];
gradK_d_0 = [0, 0, 0, 0; 
             0, 0, 0, 0];
end

function [A,gradA] = a(x) % First component for smoothing function
global par
[Psi,gradPsi_top,gradPsi,gradgradPsi_top] = psi(x);
[K_d_0, gradK_d_0] = k_d_0(x);
Lgpsi = gradPsi_top;
A = Lgpsi*K_d_0 + par.alpha*Psi;
gradA = (K_d_0.')*gradgradPsi_top + Lgpsi*gradK_d_0 + par.alpha*gradPsi;
end

function [B,gradB] = b(x) % Second component for smoothing function
[Psi,gradPsi_top,gradPsi,gradgradPsi_top] = psi(x);
Lgpsi = gradPsi;
B = norm(Lgpsi)^2;
gradB = 2*gradPsi_top*gradgradPsi_top;
end

function [HalfSontag, gradHalfSontag] = half_sontag(x) % Smoothing function and its derivative
global par
[A,gradA] = a(x);
[B,gradB] = b(x);
if B == 0
    HalfSontag = 0;
    gradHalfSontag = 0;
else
    HalfSontag = (-A+sqrt(A^2 + par.sigma*B^2))/(2*B);
    gradHalfSontag = ((-gradA+(A*gradA+par.sigma*B*gradB)/sqrt(A^2 + par.sigma*B^2))*2*B-2*gradB*(-A+sqrt(A^2 + par.sigma*B^2)))/(4*B^2);

end
end

function [K_0, gradK_0] = k_0(x) % Virtual controller kappa and its derivative
global par
[Psi,gradPsi_top,gradPsi,gradgradPsi_top] = psi(x);
[K_d_0, gradK_d_0] = k_d_0(x);
[HalfSontag, gradHalfSontag] = half_sontag(x);
Lgpsi = gradPsi_top;
K_0 = K_d_0 + HalfSontag*Lgpsi.';
gradK_0 = gradK_d_0 + (Lgpsi.')*gradHalfSontag + HalfSontag*gradgradPsi_top;
end

function [h,gradh] = ABC(x) % Activated Backstepping CBF
global par
[Psi,gradPsi_top,gradPsi,gradgradPsi_top] = psi(x);
[Top, gradTop] = top_layer(x);
[K_0, gradK_0] = k_0(x);
h = Psi - requ_mine(-(1/(2*par.mu))*gradPsi_top*(Top-K_0));
gradh = gradPsi - 2*relu_mine(-(1/(2*par.mu))*gradPsi_top*(Top-K_0))*(-1/(2*par.mu))*(gradPsi_top*(gradTop-gradK_0)+((Top-K_0)')*gradgradPsi_top);
end

function out = k_d(x) % Desired controller for the safe controller based on ABC
global par
out = [-par.Keta*x(2)-par.Ktheta*sin(x(3)); par.Kv*(par.vd-x(4))];
end

function out = u_ABC(x)
global par
[h,gradh] = ABC(x);
H = [1, 0; 0, 0.15]; % Gamma weights
F = [];
Aeq = [];
beq = [];
x0 = [0; 0];
A = -gradh*g(x);
b = par.alpha_2*h+gradh*f(x)+gradh*g(x)*k_d(x);
options = optimoptions('quadprog','Display','off','Algorithm','active-set');
out = quadprog(H,F,A,b,Aeq,beq,[],[],x0,options) + k_d(x);
end

function out = relu_mine(x) % ReLU function
if x <= 0
    out = 0;
else
    out = x;
end
end

function out = requ_mine(x) % ReQU function
if x <= 0
    out = 0;
else
    out = x^2;
end
end