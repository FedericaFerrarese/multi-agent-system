% main file 
clc
close all
clear all

% Fix parameters
N     = 500;                % number of preys
dt    = 0.1;                % time step
tf    = 10; t0 = 0;         % time interval 
steps = (tf-t0)/dt + 1;     % number of steps


par.R =0.1;                 % repulsion between predators
par.p      = 3;             % distance influence
par.a      = 1;             % attraction between preys
par.b      =1.5;            % repulsion between preys and predators 


% selection of the initial data 
disp ('[1] 1 predator')
disp('[2] 2 predators')
P = input('initial data 1-2: ');

% one predator
if P==1
% attraction between predators and preys 
%par.c    = 1; % esce
%par.c = 3; % entra 
par.c = 6; % entra e oscilla
% par.c = 10; % aumentano le oscillazioni 
%par.c=5;

% initial data
% preys positions
x           =  -0.5+rand(N,1);        
y           =  -1.5+3*rand(N,1);

% predator position
zx          = -1+2*rand(P,1);     
zy          =  -1+2*rand(P,1);     
end

% two predators
if P>1
% attraction between predators and preys 
    par.c =[4 4.1];
% preys positions    
    x = -0.5+rand(N,1);        
    y = -1.5+3*rand(N,1);
% predators positions
    zx = [-1 ; 1];
    zy = [0;0];
end

% plot of initial data
figure
plot(x,y,'bo','LineWidth',2.5)
hold on
plot(zx,zy,'ro','LineWidth',2.5)
title([' Preys-predator model for t=' num2str(1),',for a=' num2str(par.a),...
',b =' num2str(par.b), ',c =' num2str(par.c),',p =' num2str(par.p)])
axis([-1,1,-1,1])
axis('equal')
xlabel('x','FontSize',12,'FontWeight','bold')
ylabel('y','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold')


% initial data 
z0 = [x;y;zx;zy];

% set the tollerances 
opts = odeset('RelTol',1e-2,'AbsTol',1e-4);

% numerical solution 
[t,z]= ode45(@(t,z) modelAR_2D(t,z,par,P),linspace(t0,tf,steps),z0,opts);

% time loop
figure
for i = 1:steps-1
        
        % divide the vector z in each component
        x  = z(i,1:N);
        y  = z(i,N+1:2*N);
        zx = z(i,2*N+1:2*N+P);
        zy = z(i,2*N+P+1:2*N+2*P);
        
        % plot the time evolution 
        plot(x,y,'bo','LineWidth',2.5)
        hold on 
        plot(zx,zy,'ro','LineWidth',2.5)
         axis([-5,5,-5,5])
        axis('equal')
        title([' Preys-predator model for t=' num2str(i+1),',for a=' num2str(par.a),...
        ',b =' num2str(par.b), ',c =' num2str(par.c),',p =' num2str(par.p)])
        pause(0.1)
        hold off    
        drawnow
end
xlabel('x','FontSize',12,'FontWeight','bold')
ylabel('y','FontSize',12,'FontWeight','bold')
 
set(gca,'FontSize',12,'FontWeight','bold')

