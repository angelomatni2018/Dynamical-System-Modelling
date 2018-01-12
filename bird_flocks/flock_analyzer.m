function flock_analyzer

seed=17; rng(seed);

% system parameters

% parameters for all choices (can be overwritten)
parameters.dt=0.1;
parameters.tmax=100;
opt.use_noise = 0;

opt.method='rk4';

% plotting options
% plots movie only for time > plot_movie
% opt.plot_movie=50;
opt.plot_movie = inf;

% plot graphs if option set to 1
% opt.plot_graphs = 1;
opt.plot_graphs = 0;

% allow a sequence of dt for convergence test and loop over it
%dt=dt0*2.^-[0:0];

% type of parameters to use (based on part of lab)
global parameter_choice;

parameter_choice = 12;

switch parameter_choice
    case 1 % for vortices
        parameters.radius=1;
        parameters.N_boids=10;
        parameters.L=20;
        parameters.radat=2.5;
        parameters.radrep=1.4;
        parameters.delta_rep=0.45;
        parameters.noise=0;
        parameters.alpha=0.5;%0.25;
        parameters.beta=2;%2;
        init_cond_choice=1;
    case 2 % for citcular trajectory test
        parameters.radius=.5;
        parameters.N_boids=2;
        parameters.L=20;
        parameters.radat=2.5;
        parameters.radrep=1.4;
        parameters.delta_rep=0.45;
        parameters.noise=0;
        parameters.alpha=0.5;%0.25;
        parameters.beta=2;%2;
        init_cond_choice=3;
    case 3 % for convergence
        parameters.tmax=10;
        parameters.radius=1;
        parameters.N_boids=10;
        parameters.L=20;
        parameters.radat=2.5;
        parameters.radrep=1.4;
        parameters.delta_rep=0.45;
        parameters.noise=0;
        parameters.alpha=0.5;%0.25;
        parameters.beta=2;%2;
        init_cond_choice=2;
    case 10 % section 5 part a
        parameters.dt=.05;
        parameters.radius=.25;
        parameters.N_boids=30;
        parameters.L=10;
        parameters.radat=1;
        parameters.radrep=1;
        parameters.delta_rep=0.45;
        parameters.noise=0;
        parameters.alpha=0;
        parameters.beta=0;
        init_cond_choice=1;
        opt.plot_graphs = 1;
    case 11 % section 5 part b
        parameters.alpha=0;
        parameters.beta=0;
        parameters.radat=2.5;
        parameters.radrep=1.4;
        parameters.delta_rep=0.45;
        parameters.N_boids=320;
        parameters.L=30;
        parameters.radius=2.5;
        opt.use_noise = 1;
        init_cond_choice=1;
   case 12 % section 5 part c
        parameters.L=50;
        parameters.N_boids=100;
        parameters.noise=0;
        parameters.alpha=.5;
        parameters.beta=2;
        parameters.radat=2.5;
        parameters.radrep=1.4;
        parameters.delta_rep=0.45;
        parameters.dt=.05;
        parameters.tmax=100;
        parameters.radius=1;
        init_cond_choice=1;
        opt.plot_graphs = 1;
end

switch init_cond_choice
    case 1 % random across whole system
        init_cond.x0=parameters.L*rand(parameters.N_boids,1);
        init_cond.y0=parameters.L*rand(parameters.N_boids,1);
        init_cond.vx0=2*(rand(parameters.N_boids,1)-.5);
        init_cond.vy0=2*(rand(parameters.N_boids,1)-.5);
        % Normalize
        vnorm=sqrt(init_cond.vx0.^2+init_cond.vy0.^2);
        init_cond.vx0=init_cond.vx0./vnorm;
        init_cond.vy0=init_cond.vy0./vnorm;
    case 2 % convergence square
        sq_x = parameters.L*rand(parameters.N_boids,1);
        sq_x = 7/16*parameters.L*ones(size(sq_x)) + 1/8*parameters.L*sq_x;
        sq_y = parameters.L*rand(parameters.N_boids,1);
        sq_y = 7/16*parameters.L*ones(size(sq_y)) + 1/8*parameters.L*sq_y;
        init_cond.x0=sq_x;
        init_cond.y0=sq_y;
        init_cond.vx0=2*(rand(parameters.N_boids,1)-.5);
        init_cond.vy0=2*(rand(parameters.N_boids,1)-.5);
        % Normalize
        vnorm=sqrt(init_cond.vx0.^2+init_cond.vy0.^2);
        init_cond.vx0=init_cond.vx0./vnorm;
        init_cond.vy0=init_cond.vy0./vnorm; 
    case 3 % for circular trajectory test
        init_cond.x0=[10,12];
        init_cond.y0=[10,10];
        init_cond.vx0=[0,0];
        sq_dist=4;
        a_term = parameters.alpha*exp(-sq_dist/parameters.radat^2);
        b_term = parameters.beta*exp(-sq_dist/parameters.radrep^2)/...
                    (sq_dist+parameters.delta_rep^2);
        force_att_rep= a_term - b_term;
        init_cond.vy0=[0,2*sqrt(force_att_rep)];
end

switch parameter_choice
    case 3
        num_runs = 10;
        finals = NaN(num_runs,4*parameters.N_boids);
        dts = NaN(num_runs + 1);
        dts(1) = .25;
        for ii = 1:size(finals,1)
            dts(ii+1) = dts(ii) / 2;
            parameters.dt = dts(ii);
            [final,~] = flock(parameters,init_cond,opt);
            finals(ii,:) = final(end,:);
        end
        errors = diff(finals,1);
        error_sums = zeros(num_runs - 1);
        n_boids = parameters.N_boids;
        for ii = 1:n_boids
            error_sums = error_sums + errors(:,ii).^2 + ...
                errors(:,n_boids+ii).^2 + errors(:,2*n_boids+ii).^2 + ...
                errors(:,3*n_boids+ii).^2;
        end
        error_sums = sqrt(error_sums) / 4 / n_boids;
        loglog(dts(1:size(error_sums,1)),error_sums,...
            '-ro','LineWidth',2,...
            'MarkerSize',14,'MarkerFaceColor',[1,0,0]);
    case 11 % Section 5 part b
        parameters.dt=.5;
        parameters.tmax=100;
        parameters.noise=0;
        noises = linspace(0,.6,10);
        order_params = NaN(size(parameters.noise));
        for ii = 1:size(noises,2)
            disp('beginning run')
            parameters.noise = noises(ii);
            [~,order_params(ii)]=flock(parameters,init_cond,opt);
        end
        plot(noises, order_params);
    otherwise
        [xfinal,order_parameter]=flock(parameters,init_cond,opt);
end

end % flock_wrapper

%%%%%%%%%%%%%%%%%%%%%%%

function [xfinal,order_parameter]=flock(parameters,init_cond,opt)

% the boids fly with constant speed and interact with each other with a
% strength that decays with distance like a Gaussian with width radius.
% Their velocity tends towards the average velocity of the neighboring
% boids (weighted by distance with the Gaussian mentioned above). In
% addition for non-zero alpha they attract each other with a force
% decaying like a Gaussian with width radat and
% for non-zero beta they repel each other again with a Gaussian with width
% radrep
tic
global alpha beta radat radrep delta_rep N_boids L radius

dt=parameters.dt;
tmax=parameters.tmax;
radius=parameters.radius;
N_boids=parameters.N_boids;
L=parameters.L;
radat=parameters.radat;
radrep=parameters.radrep;
delta_rep=parameters.delta_rep;
noise=parameters.noise;
alpha=parameters.alpha;
beta=parameters.beta;
x0=init_cond.x0;
y0=init_cond.y0;
vx0=init_cond.vx0;
vy0=init_cond.vy0;
%%%%%%%%%%%%%%%%%%%

% preallocate the memory for the positions to speed up the code
numsteps=tmax/dt+1;
angle=NaN(numsteps,N_boids);
Y=NaN(numsteps,4*N_boids);
V=NaN(numsteps);
time=NaN(numsteps,1);

% indices for positions and velocities
x_index=1:N_boids;
y_index=N_boids+1:2*N_boids;
vx_index=2*N_boids+1:3*N_boids;
vy_index=3*N_boids+1:4*N_boids;

xp=x0;
yp=y0;
if (opt.plot_movie>0 && opt.plot_movie ~= inf) % set-up for movie
    % Using preallocated angle var:
    % [~,HandleM,HandleF]=set_up_movie(xp,yp,angle(1,:));
    % Not using preallocated angle var:
    [~,HandleM,HandleF]=set_up_movie(xp,yp,NaN(1,N_boids));
end

% start integration
Y(1,x_index)=x0;
Y(1,y_index)=y0;
Y(1,vx_index)=vx0;
Y(1,vy_index)=vy0;
time(1)=0;

% integration loop
k=1;
while (time(k) < tmax)
    if (time(k)+dt>tmax) % make sure that last time step hits tmax exactly
        dt=tmax-time(k);
    end
    switch opt.method
        case 'euler' % forward Euler
            [time(k+1),Y(k+1,:)]=euler(time(k),Y(k,:),dt);
        case 'rk4'  % RK 4
            [time(k+1),Y(k+1,:)]=rk4(time(k),Y(k,:),dt);
    end
    
    % add a bit of noise to velocities
    if (opt.use_noise)
        % determine order parameter
        order_x = 1/N_boids * sum(Y(k,vx_index));
        order_y = 1/N_boids * sum(Y(k,vy_index));
        V(k) = sqrt(order_x^2 + order_y^2);
        % add noise to velocity
        Y(k+1,vx_index) = Y(k+1,vx_index) + sqrt(dt)*noise*randn(size(vx_index));
        Y(k+1,vy_index) = Y(k+1,vy_index) + sqrt(dt)*noise*randn(size(vy_index));
        % normalize
        magnitudes = sqrt(Y(k+1,vx_index).^2 + Y(k+1,vy_index).^2);
        Y(k+1,vx_index) = Y(k+1,vx_index) ./ max(magnitudes,1e-10);
        Y(k+1,vy_index) = Y(k+1,vy_index) ./ max(magnitudes,1e-10);
    end
    
    k=k+1;
    % periodic boundary conditions
    xp=mod(Y(k,x_index),L);
    yp=mod(Y(k,y_index),L);
    angle(k,:)=360*atan2(Y(k,vy_index),Y(k,vx_index))/(2*pi);
    if (time(k)>opt.plot_movie) % update plot
        set(HandleM,'XData',xp,'MarkerFaceColor','r');
        set(HandleM,'YData',yp);
        set(HandleF,'CData',angle(k,:)');
        drawnow
    end
end
kmax=k;

toc

order_parameter = NaN;
% need to computer order parameter
if (opt.use_noise)
    V_avg=sum(V(round(kmax/2):kmax-1))/length(round(kmax/2):kmax-1);
    order_parameter = V_avg;
end

% output trajectories

if (opt.plot_graphs == 1)
    figure(2);
    plot(Y(:,x_index),Y(:,y_index),'-');
    title(' Trajectories of each boid');
    xlabel(' x-position');
    ylabel(' y-position');

    figure(3);
    plot(time,angle,'o');
    title(' Angles of trajectories of each boid');
    xlabel(' Time (in seconds)');
    ylabel(' Angle (in degrees)');
end

xfinal=Y(kmax,:);

end % flock


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dY] = F(time,Y)

global alpha beta radat radrep delta_rep N_boids L radius parameter_choice

if parameter_choice == 2
    dY = F_traj(time,Y);
    return;
end

Fvx=NaN(1,N_boids);
Fvy=NaN(1,N_boids);

xp=mod(Y(1:N_boids),L);
yp=mod(Y(N_boids+1:2*N_boids),L);
vx=Y(2*N_boids+1:3*N_boids);
vy=Y(3*N_boids+1:4*N_boids);

for i_boid=1:N_boids
    % distances of all boids from boid i_boid
    distancex=xp(i_boid)*ones(1,N_boids)-xp;
    distancey=yp(i_boid)*ones(1,N_boids)-yp;
    % periodic boundary conditions
    distancex(distancex>L/2)=distancex(distancex>L/2)-L;
    distancex(distancex<-L/2)=distancex(distancex<-L/2)+L;
    distancey(distancey>L/2)=distancey(distancey>L/2)-L;
    distancey(distancey<-L/2)=distancey(distancey<-L/2)+L;
    distancesq=max(distancex.^2+distancey.^2,.1e-10);  % avoid dividing by 0 in the force term below
    weights=exp(-distancesq/radius^2);
    %weighted average of the velocities
    vxaverage_weighted=vx*weights'/sum(weights);
    vyaverage_weighted=vy*weights'/sum(weights);
    vavnorm=sqrt(vxaverage_weighted^2+vyaverage_weighted^2);
    force_att_rep=alpha*exp(-distancesq'/radat^2)-...
        beta*exp(-distancesq'/radrep^2)./(distancesq'+delta_rep^2*ones(N_boids,1));
    Fvx(1,i_boid)=vxaverage_weighted/vavnorm-vx(i_boid)-distancex*force_att_rep;
    Fvy(1,i_boid)=vyaverage_weighted/vavnorm-vy(i_boid)-distancey*force_att_rep;
end

dY=[vx';vy';Fvx';Fvy'];

end % F

function [dY] = F_traj(time,Y)

global alpha beta radat radrep delta_rep N_boids L

xp=mod(Y(1:N_boids),L);
yp=mod(Y(N_boids+1:2*N_boids),L);
vx=Y(2*N_boids+1:3*N_boids);
vy=Y(3*N_boids+1:4*N_boids);

Fvx=NaN(1,N_boids);
Fvy=NaN(1,N_boids);

for i_boid=1:N_boids
    % distances of all boids from boid i_boid
    distancex=xp(i_boid)*ones(1,N_boids)-xp;
    distancey=yp(i_boid)*ones(1,N_boids)-yp;
    % periodic boundary conditions
    distancex(distancex>L/2)=distancex(distancex>L/2)-L;
    distancex(distancex<-L/2)=distancex(distancex<-L/2)+L;
    distancey(distancey>L/2)=distancey(distancey>L/2)-L;
    distancey(distancey<-L/2)=distancey(distancey<-L/2)+L;
    distancesq=max(distancex.^2+distancey.^2,.1e-10);  % avoid dividing by 0 in the force term below
    force_att_rep=alpha*exp(-distancesq'/radat^2)-...
        beta*exp(-distancesq'/radrep^2)./(distancesq'+delta_rep^2*ones(N_boids,1));
    Fvx(1,i_boid)=-distancex*force_att_rep;
    Fvy(1,i_boid)=-distancey*force_att_rep;
end

% Force the first boid to be stationary
Fvx(1,1) = 0;
Fvy(1,1) = 0;

dY=[vx';vy';Fvx';Fvy'];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [time,Y]=euler(time,Y,dt)
dY=F(time,Y);
Y=Y+dt*dY';
time=time+dt;
end

function [time,Y]=rk4(time,Y,dt)
k1 = F(time,Y);
k2 = F(time + .5*dt, Y + .5*dt*k1');
k3 = F(time + .5*dt, Y + .5*dt*k2');
Y = Y + dt/6 * (k1' + 2*k2' + 2*k3' + F(time + dt, Y + dt*k3')');
time = time + dt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [HandleFig,HandleM,HandleF]=set_up_movie(xp,yp,angle)
global L N_boids
HandleFig=figure(1);
screen=get(0,'ScreenSize');
width=.70*screen(3);
height=.60*screen(4);
set(HandleFig,'Position',[.15*screen(3),.20*screen(4),width,height]);
subplot(1,6,[1:4])
HandleM=plot(xp,yp,'o','MarkerSize',5);
title(' Flocks of Boids')
set(gcf,'DoubleBuffer','on');
axis([0 L 0 L]);
axis square
xlabel(' x-Position of Boids')
ylabel(' y-Position of Boids')
%figure(10)
subplot(1,6,[5:6])
ylabel(' Number of Boid')
HandleF=imagesc(angle');
caxis([-180 180])
axis([0.5 1.5 0 N_boids])
colorbar
end
