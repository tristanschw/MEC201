%% Project 4: Genetic Algorithm Optimization
%% Genetic Algorithm
children = 6;
parents = 6;
S = 20;
G = 100;
theta_desired = 650;
TOL = 1*10^(-6);
w1 = 1000;
w2 = 100; 

% Definitions
Pi_a(theta)=w1*((theta_max - theta_des)/theta_des)^2 + w2*d_f; 


while (i <= S)
    D = (0.5 + 0.4 .* rand);
    J = 10^7 .* rand;
    Lambda = Lambda(index);
    Lambda_2 = Lambda_2(index);
    [max_theta,d_f] = densification(D,J);
    pi_val = densification(D,J);
    new_pi = [new_pi pi_val];
    i = i + 1; 
end
[new_pi, ind] = sort(new_pi); %Sort new pi
PI(1,:) = new_pi;             %Store new pi
PI_min = [PI_min min(new_pi)];%Store cost of best performer                    
PI_avg = [PI_avg mean(new_pi)];%Store average cost of performers
Lambda = Lambda(ind);          %Store Lambda
Orig(1,:) = ind;               %Store indices of the first generation

for s=1:S
    new_pi(1,s) =densification(Lambda(s,1));
end

while PI_min(end)>TOL_GA && g<=G %While the generation limit has not been exceeded and the cost is greater than TOL
     g=g+1; 
     for u = 1:2:P %For Populating 
         
         % Assign random numbers to Phi
         Phi_1 = randi([0 100]);
         Phi_2 = randi([-100 0]);
         
         %Generating New Children using Equation 5
         D_c1 = Phi_1*new_pi(u) + (1-Phi_1)*new_pi(u+1); 
         J_c2 = Phi_2*new_pi(u) + (1-Phi_2)*new_pi(u+1);
         
         Lambda(P+u,dv) = D_c1;
         Lambda(P+u+1,dv) = J_c2;
         pi_c1 = functio(D_c1);
         pi_c2 = functio(J_c2);
         new_pi(dv,P+u) = pi_c1;
         new_pi(dv,P+u+1) = pi_c2;
     end
     for n=S-2*P:S
         Lambda = -20+(20+20)*rand(S,1);
         new_pi(1,n) = densification(Lambda(n,1));
     end
     [new_pi, ind] = sort(new_pi); %Sort new pi
     PI(1,:) = new_pi;             %Store new pi
     PI_min = [PI_min min(new_pi)];%Store cost of best performer
     PI_avg = [PI_avg mean(new_pi)];%Store average cost of performers
     Lambda = Lambda(ind);          %Store Lambda
     Orig(1,:) = ind;
end
