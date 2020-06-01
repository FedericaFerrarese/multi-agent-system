% % main file 
clc
close all
clear all
% Fix parameters
N     =50;                % number of preys 
dt    = 0.1;               % time step
tf    = 10; t0 = 0;        % time interval 
steps = (tf-t0)/dt + 1;    % number of steps

par.p      = 3;             % distance influence
par.a      = 1;             % attraction between preys
par.b      =0.2;            % repulsion between preys and predator 
par.c    = 1.8;             % attraction between predator and preys



% initial data
% preys positions
x           =  -0.5+rand(N,1);        
y           =  -1.5+3*rand(N,1);

% predator position
zx          = -1+2*rand;     
zy          =  -1+2*rand;  

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
z=z0;



Rz   = sqrt((x-zx').^2 + (y-zy').^2);  % matrix of distances
N1(1)=N; % number of preys at initial time
t(1)=0;  % initial time

B   = @(r) par.b./(r.^2);        % preys predators repulsion
A  = @(r) r.^(-2)-par.a;         % preys preys repulsion and attraction
C= @(r)r.^(-par.p);              % attraction between predators and preys


figure
for i1= 1:steps-1

        V=Rz>0.5; % truth vector
        % survivor preys
        x=x.*V;    
        y=y.*V;
        % remove the death preys 
        x(x==0)=[];
        y(y==0)=[];
        
        N=sum(V); % number of survivors 
        N1(i1+1)=N; % number of preys at time step i1+1 
        t(i1+1)=i1; % time 

        
        Rx   = sqrt((x-x').^2 + (y-y').^2+eye(N));   % matrix of distances
        Rz   = sqrt((x-zx').^2 + (y-zy').^2);        % matrix of distances


        FA   = A(Rx);  % repulsion and attraction force between preys
        FB = B(Rz);    % repulsion force between preys and predators 
        FC = C(Rz);    % attraction force between predators and preys

        % preys force 
        vx = x; vy = y;
        for i=1:N
            vx (i) =  ((sum(FA(i,:)'.*(x(i)-x)))./N + FB(i).*(x(i)-zx));
            vy (i) =  (sum(FA(i,:)'.*(y(i)-y)))./N + FB(i).*(y(i)-zy);
        end
        
        % predator force 
        vzx = sum(par.c*FC.*(x-zx))/N;
        vzy = sum(par.c*FC.*(y-zy))/N;
        
        
        % Euler method 
        x = x+ dt*vx;
        y = y+dt*vy;
        zx=zx+dt*vzx;
        zy=zy+dt*vzy;
        
        % time evolution 
        plot(x,y,'bo','LineWidth',2.5)
        hold on 
        plot(zx,zy,'ro','LineWidth',2.5)
        axis([-5,5,-5,5])
        title([' Preys-predator model for t=' num2str(i1+1),',for a=' num2str(par.a),...
        ',b =' num2str(par.b), ',c =' num2str(par.c),',p =' num2str(par.p)])
        pause(0.1)
        hold off    
        drawnow
 end
xlabel('x','FontSize',12,'FontWeight','bold')
ylabel('y','FontSize',12,'FontWeight','bold') 
set(gca,'FontSize',12,'FontWeight','bold')


% preys death 
figure
plot(t,N1,'b','linewidth',3)
title('Preys decay in time')
xlabel('t','FontSize',12,'FontWeight','bold')
ylabel('N','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold')
axis([0 50 0 50])