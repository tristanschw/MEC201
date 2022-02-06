% Tristan Schwab, UC BERKELEY ME, tristanschwab@berkeley.edu
% Prepared for ME201, Project One

save_images = 1; % set this to 1 to save outputs
output_folder = "./output_images";
if save_images &&  ~(isfolder(output_folder))
    mkdir(output_folder);
end

%% Problem One, Graphing Objective Functions

% House-Cleaning
close all % close all open figures
clc % clear command window
clear

% Definitions
x = linspace(-20,20);
P1=x.^2;
P2=(x+(pi/2)*sin(x)).^2;
g=diff(P1);
h=diff(P2);

%Plotting
figure("name", "Objective Functions")  % Generate new named figure for plot
plot(x, P1,'b', 'LineWidth', 3) % use 'LineWidth' attribute to increase line
% thickness. Makes plot lines a lot easier to see
hold on
plot(x,P2,'g', 'LineWidth',3)
ylabel('$\Pi_{A,B}$','Interpreter','latex')
xlabel("x")
title('Original Objective Functions') % Always include title, axis labels with 
set(gca,'FontSize', 26); % Use this command to increase text size on plot
axis tight % make axis crop out any dead space
x0 = 0; y0 = 0; width = 800; height = 500;
set(gcf,'units','pixels','position',[x0,y0,width,height]) 

%% Algorithm Two for myNewton

%House-Cleaning
close all % close all open figures
clc % clear command window
clear

%Definitions
syms x Pi_a(x) Pi_b(x)
Pi_a(x)=x.^2; Pi_b(x)=(x+(pi/2)*sin(x)).^2;
f_a = diff(Pi_a,x); f_b = diff(Pi_b,x);
df_a = matlabFunction(hessian(Pi_a,x)); 
df_b = matlabFunction(hessian(Pi_b,x));

%Initial Conditions
TOL = 10^(-8);
maxit = 20;
g=[-1,0,1];
hist_a1=zeros(1,length(g));
hist_b1=zeros(1,length(g));
for k = g
    x0 = 2*10^(k);
    %Evaluate Pi_A
    [sol_a, its_a, hist_a] = myNewton(matlabFunction(f_a), df_a, x0, TOL, maxit);
    %Pi_a=[Pi_a,sol_a];

    %Evaluate Pi_B
    [sol_b, its_b, hist_b] = myNewton(matlabFunction(f_b), df_b, x0, TOL, maxit);
    %Pi_b=[Pi_b,sol_b];
    
    % Plotting 
    figure("name","Pi A");
    plot(1:1:its_a+1, subs(Pi_a,hist_a),'g','LineWidth',3);
    title(['\Pi_A vs Iterations Using X0=' num2str(x0)])    
    xlabel('Iterations')
    ylabel('Pi_a')
    hold on 
    figure("name","Pi B");
    plot(1:1:its_b+1, subs(Pi_b,hist_b),'b','LineWidth',3);
    title(['\Pi_B vs Iterations Using X0=' num2str(x0)])  
    xlabel('Iterations')
    ylabel('Pi_b')
end

%% Genetic Algorithm

% define variables

S = 50;
G = 100;
P = 12;
dv=1;
TOL_GA = 10^(-6);
PI = zeros(G,S);
Orig = zeros(G,S);
PI_min = [];
PI_avg = [];
new_pi = zeros(1,S);
g=1;
functio = Pi_a;
Lambda = -20+(20+20)*rand(S,1);

for s=1:S
    new_pi(1,s) =functio(Lambda(s,1));
end

[new_pi, ind] = sort(new_pi); %Sort new pi
PI(1,:) = new_pi;             %Store new pi
PI_min = [PI_min min(new_pi)];%Store cost of best performer                    
PI_avg = [PI_avg mean(new_pi)];%Store average cost of performers
Lambda = Lambda(ind);          %Store Lambda
Orig(1,:) = ind;               %Store indices of the first generation

while PI_min(end)>TOL_GA && g<=G %While the generation limit has not been exceeded and the cost is greater than TOL
     g=g+1;
     
     for u = 1:2:P %For Populating 
         
         % Assign random numbers to Phi
         Phi_1 = randi([0 100]);
         Phi_2 = randi([-100 0]);
         
         %Generating New Children using Equation 5
         x_c1 = Phi_1*new_pi(u) + (1-Phi_1)*new_pi(u+1); 
         x_c2 = Phi_2*new_pi(u) + (1-Phi_2)*new_pi(u+1);
         
         Lambda(P+u,dv) = x_c1;
         Lambda(P+u+1,dv) = x_c2;
         pi_c1 = functio(x_c1);
         pi_c2 = functio(x_c2);
         new_pi(dv,P+u) = pi_c1;
         new_pi(dv,P+u+1) = pi_c2;
     end
     for n=S-2*P:S
         Lambda = -20+(20+20)*rand(S,1);
         new_pi(1,n) = functio(Lambda(n,1));
     end
     [new_pi, ind] = sort(new_pi); %Sort new pi
     PI(1,:) = new_pi;             %Store new pi
     PI_min = [PI_min min(new_pi)];%Store cost of best performer
     PI_avg = [PI_avg mean(new_pi)];%Store average cost of performers
     Lambda = Lambda(ind);          %Store Lambda
     Orig(1,:) = ind;
end

figure("name","PI_min vs. PI_avg")
semilogy(1:g,PI_min, 1:g, PI_avg)
title('Log Plot of Average and Minimum Cost vs. Number of Generations')
xlabel('Generations')
ylabel('Cost')
legend('PI min','PI avg')


