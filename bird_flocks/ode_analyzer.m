function ode_analyzer
% solves ode with Euler or RK4 with
% time-step control and Richardson extrapolation

set(0,'DefaultAxesFontsize',20,'defaultaxeslinewidth',1.5,...
'defaultlinelinewidth',2.,'defaultpatchlinewidth',1.5,'DefaultTextFontSize',20)

clear
format long
cputime1=cputime;

global ifunction y0 adapt time_pert perturb_factor perturb_duration

marker=make_markers;

y0=1;

% Euler does not work well with adaptive dt: errors that are made early at
% small t - when the adaptive code chooses large time steps - lead to much
% larger errors later than errors that are made closer to the singularity
% that's because an early shift in y leads to a much larger shift in the
% time of the singularity than the same error in y made later (since the
% slope of the solution increases so much). But why does this not carry
% over directly to RK2 and RK4? Is the order the relevant quantity? 
% using also Richardson extrapolation in the FE adaptive makes it faster
% than the non-adaptive code (at least when tmax is not too close to the
% singularity.
% for tmax=0.995 Euler without extrapol goes for jmax=13 to error=3e-4
% (with 300,000 steps)
% for tmax=0.995 Euler with extrapol goes for jmax=13 to error=2e-8
% (with 300,000 steps)
% for tmax=0.995 RK4  without extrapol goes for jmax=9 to error=1e-10
% (with 800 steps)
% for tmax=0.995 RK4  with extrapol goes for jmax=9 to error=1e-12
% (with 800 steps)

ifunction=1;

% method: 1 = FE 2 = RK4  3 = RK2 4 = ODE45
method=2;

if (method==1)
    % FE
    p=2;
    dtfactor=10^(1/4);
    jmax=4; tol0=.1; richardson=1; adapt=1; dt_init=1;    
    %jmax=18; tol0=.1; richardson=1; adapt=0; dt_init=.1;
end
if (method==2)
    % RK4
    p=5;
    %for adaptive
    dtfactor=10^(1/1);
%    jmax=17; tol0=.01; richardson=0; adapt=1; dt_init=1;    
    jmax=8; tol0=.01; richardson=1; adapt=1; dt_init=1;    
    %jmax=2; tol0=.001; richardson=0; adapt=1; dt_init=1;
    %for perturbation
    %jmax=1; tol0=.00000001; richardson=1; adapt=1; dt_init=1;
    %jmax=1; tol0=.01; richardson=1; adapt=1; dt_init=1;
    %for non-adaptive 
    %dtfactor=10^(1/2);
    %jmax=10; tol0=.01; richardson=1; adapt=0; dt_init=.2;
end
if (method==3)
    % RK2
    p=3;
    dtfactor=2;
    jmax=10; tol0=.1; richardson=1; adapt=1; dt_init=1;
    jmax=10; tol0=.1; richardson=1; adapt=0; dt_init=.2;
end
if (method==4)
    jmax=4; tol0=.01;
end

% change tolerance briefly
time_pert=Inf;
perturb_factor=80;
perturb_duration=0.2;


% allocate memory to the arrays to speed up the code
ntmax=2000000;
% sol stores position of particle for all times
sol=NaN(ntmax,1);
% for variable time step need to keep track of deltat and t separately
deltat=NaN(ntmax,1);
time=NaN(ntmax,1);
error_cum=NaN(ntmax,1);
error_est=NaN(ntmax,1);

tmax=.9999;%.9999;
dtminimal=.1e-11;

for j=1:jmax
    tolj(j)=tol0/(dtfactor^((2-1)*(j-1)));
    dtj(j)=dt_init/(dtfactor^((2-1)*(j-1)));
end

for j=1:jmax
    tol=tolj(j);
    if (adapt==1)
        dt=dt_init;
    else
        dt=dtj(j);
    end
    time_i=0;
    nt=1;
    sol(1)=y0;
    y=y0;
    time(1)=time_i;
    deltat(1)=dt;
    error_cum(1)=0;
    error_est(1)=0;
    
    tic
    if (method==4)
        options=odeset('AbsTol',tol,'RelTol',3e-14);
        [time,sol]=ode45(@F,[0,tmax],y0,options);
        ntj(j)=length(time);
        
    else
        
        while time_i < tmax
            
            if (time_i+dt > tmax-dtminimal)
                dt=tmax-time_i;
            end
            
            if (time_i>time_pert && time_i<time_pert+perturb_duration)
                tol_pert=tol*perturb_factor;
            else
                tol_pert=tol;
            end
            
            if (method==1)
                %  [time_i,y,dt,success]=stepeuleroc(time_i,y,dt,tol);
                %  [time_i,y,dt,success]=stepeuler(time_i,y,dt,tol);
                [time_i,y,dt,success,errorestimate]=stepeuleropt(time_i,y,dt,tol_pert,tmax,richardson,adapt);
                %   [time_i,y,dt,success]=stepeuleropt(time_i,y,dt,tol,tmax);
            end
            
            if (method==2)
                [time_i,y,dt,success,errorestimate]=steprk4opt(time_i,y,dt,tol_pert,tmax,richardson,adapt);
                %  [time_i,y,dt,success]=steprk(time_i_i,y,dt,tol);
                %  [time_i,y,dt,success]=steprkoc(time_i,y,dt,tol);
                %  [time_i,y,dt,success]=steprkocwe(time_i,y,dt,tol);
            end
            
            if (method==3)
                [time_i,y,dt,success,errorestimate]=steprk2opt(time_i,y,dt,tol_pert,tmax,richardson,adapt);
                %  [time_i,y,dt,success]=steprk(time_i_i,y,dt,tol);
                %  [time_i,y,dt,success]=steprkoc(time_i,y,dt,tol);
                %  [time_i,y,dt,success]=steprkocwe(time_i,y,dt,tol);
            end

            if (success > 0)
                nt=nt+1;
                sol(nt)=y;
                time(nt)=time_i;
                deltat(nt)=dt;
                error_cum(nt)=error_cum(nt-1)+errorestimate;
                error_est(nt)=errorestimate;
            else
                if (dt<dtminimal)
                    fprintf(' dt is too small %g \n',dt)
                    break
                end
            end
            
            if (mod(nt,1000000)==0)
                fprintf(' %g ',time(nt))
            end
            
        end
        
        ntj(j)=nt;
        cput(j)=toc; 
        % these assingments would resize the variables and I assume that mnemory wouldhave to be realocated if in the next run the nt is larger
        %time=time(1:nt);
        %sol=sol(1:nt);
        %deltat=deltat(1:nt);
        %error_cum=error_cum(1:nt);        
        %error_est=error_est(1:nt);
        
    end
    uexact_t=yexact(time(1:nt),y0);
    error_t=abs(uexact_t(1:nt)-sol(1:nt));
    if (adapt==1)
        error_local(2:nt)=yexact(time(2:nt)-time(1:nt-1),sol(1:nt-1))-sol(2:nt);
        figure(20)
%        loglog(1-time(2:nt),abs(error_local(2:nt)),marker{j},1-time(2:nt),abs(error_est(2:nt)/2^p),marker{100+j})
% the estimated error for dt/2 seems incorrect, should be divided by 2^(p-1)
%        loglog(1-time(2:nt),abs(error_local(2:nt)),marker{j},1-time(2:nt),abs(error_est(2:nt)/2^(p-1)),marker{100+j})
        loglog(1-time(2:nt),abs(error_local(2:nt)')./deltat(2:nt),marker{j},1-time(2:nt),abs(error_est(2:nt)/2^(p-1))./deltat(2:nt),marker{100+j})
        xlabel(' 1 - Time')
        ylabel(' Error')
        legend(' True Local Error/dt for 1 step',' Estimated Local Error/dt');
        title([' Richardson = ',num2str(richardson)])
        hold all
        figure(21)
%        loglog(1-time(2:nt),abs(error_local(2:nt)),marker{j},1-time(2:nt),abs(error_est(2:nt)/2^p),marker{100+j})
% the estimated error for dt/2 seems incorrect, should be divided by 2^(p-1)
%        loglog(1-time(2:nt),abs(error_local(2:nt)),marker{j},1-time(2:nt),abs(error_est(2:nt)/2^(p-1)),marker{100+j})
        loglog(1-time(2:nt),abs(error_local(2:nt)'),marker{j},1-time(2:nt),abs(error_est(2:nt)/2^(p-1)),marker{100+j})
        xlabel(' 1 - Time')
        ylabel(' Error')
        legend(' True Local Error for 1 step',' Estimated Local Error');
        title([' Richardson = ',num2str(richardson)])
        hold all
    end
    if (adapt==0)
        figure(10)
        plot(time(1:nt),sol(1:nt),'o')
        xlabel(' Time'); ylabel(' y');
        hold all
    end
    
    if (method<3)
        figure(11)
        %        semilogy(time(1:nt),max(abs(error_cum(1:nt)),1e-10),'-o',time(1:nt),error_t(1:nt),'-x');
        %        semilogy(time(1:nt),abs(error_cum(1:nt)),'-o',time(1:nt),error_t(1:nt),'-x');
        loglog(1-time(1:nt),error_t(1:nt),marker{j});
        xlabel(' 1 - Time','FontSize',20)
        ylabel('Factual Total Error','FontSize',20)
        hold all
        if (adapt==1)
            figure(12)
            %        semilogy(time(1:nt),max(abs(error_cum(1:nt)),1e-10),'-o',time(1:nt),error_t(1:nt),'-x');
            %        semilogy(time(1:nt),abs(error_cum(1:nt)),'-o',time(1:nt),error_t(1:nt),'-x');
            loglog(1-time(1:nt),abs(error_cum(1:nt)),marker{j});
            xlabel(' 1 - Time')
            ylabel('Cumulative Estimated Error')
            hold all
            figure(13)
            %        semilogy(time(2:nt),max(abs(error_est(2:nt)),1e-10)./deltat(2:nt),'-x');
            %        loglog(1-time(2:nt),max(abs(error_est(2:nt)),1e-10)./deltat(2:nt),'-x');
            loglog(1-time(2:nt),abs(error_est(2:nt))./deltat(2:nt),'-x');
            xlabel(' 1 - Time')
            ylabel(' Estimated Error/dt')
            hold all
            figure(17)
            %        semilogy(time(2:nt),max(abs(error_est(2:nt)),1e-10)./deltat(2:nt),'-x');
            %        loglog(1-time(2:nt),max(abs(error_est(2:nt)),1e-10)./deltat(2:nt),'-x');
            size(time(2:nt))
            size(cumsum(error_local))
           
            loglog(1-time(1:nt),cumsum(error_local)',marker{j});
            xlabel(' 1 - Time')
            ylabel(' Cumulative Local Error')
            hold all
            %figure(14)
            %        semilogy(time(2:nt),max(abs(error_est(2:nt)),1e-10),'-x');
            %        loglog(1-time(2:nt),max(abs(error_est(2:nt)),1e-10),'-x');
            %loglog(1-time(2:nt),abs(error_est(2:nt)),'-x');
            %xlabel(' 1-Time')
            %ylabel(' Estimated Error')
            %hold all
        end
    end
        %figure(15)
        %plot(time(1:nt),sol(1:nt),'-o',time(1:nt),uexact_t(1:nt),'-')
        %hold all
    %end
    
    ufinal=sol(ntj(j));
    uexact=yexact(time(ntj(j)),y0);
    numerror(j)=abs(uexact-ufinal);
    if (adapt==1)
        fprintf(' T = %g uf = %g uex = %g  err = %g run = %g tol = %g nt = %g ctime = %g \n',time(nt),ufinal,uexact,numerror(j),j,tolj(j),ntj(j),cput(j));
    else
        fprintf(' T = %g uf = %g uex = %g  err = %g run = %g dt = %g nt = %g ctime = %g \n',time(nt),ufinal,uexact,numerror(j),j,dtj(j),ntj(j),cput(j));
    end
   
    deltat2=deltat(floor(ntj(j)/4):ntj(j)-1);
    dtmax(j)=max(deltat2);
    dtmin(j)=min(deltat2);
         
    %figure(1);
    %plot(t(1:nt),sol(1:nt));
    %ylabel('u');
    %figure(2);
    %plot(t(1:nt),deltat(1:nt),'+');
    %ylabel('Delta t');
    %xlabel('t')
    if (adapt==1)
        figure(4);
        semilogy(time(1:ntj(j)),deltat(1:ntj(j)),'+');
        ylabel('Delta t');
        xlabel('t')
        hold all
    end

end


% now compare with one run with fixed time step
time_i=0;
dt_fixed=tmax/ntj(jmax);
nt=1;
sol(1)=y0;
y=y0;
time(1)=time_i;
tic
while time_i < tmax

    if (time_i+dt_fixed > tmax-dtminimal)
        dt_fixed=tmax-time_i;
    end
    if (method==1)
        %  [time_i,y,dt,success]=stepeuleroc(time_i,y,dt,tol);
        %  [time_i,y,dt,success]=stepeuler(time_i,y,dt,tol);
        [time_i,y,dt,success,errorestimate]=stepeuleropt(time_i,y,dt_fixed,tol,tmax,richardson,0);
        %   [time_i,y,dt,success]=stepeuleropt(time_i,y,dt,tol,tmax);
    end

    if (method==2)
        [time_i,y,dt,success,errorestimate]=steprk4opt(time_i,y,dt_fixed,tol,tmax,richardson,0);
        %  [time_i,y,dt,success]=steprk(time_i_i,y,dt,tol);
        %  [time_i,y,dt,success]=steprkoc(time_i,y,dt,tol);
        %  [time_i,y,dt,success]=steprkocwe(time_i,y,dt,tol);
    end

    if (method==3)
        [time_i,y,dt,success,errorestimate]=steprk2opt(time_i,y,dt_fixed,tol,tmax,richardson,0);
        %  [time_i,y,dt,success]=steprk(time_i_i,y,dt,tol);
        %  [time_i,y,dt,success]=steprkoc(time_i,y,dt,tol);
        %  [time_i,y,dt,success]=steprkocwe(time_i,y,dt,tol);
    end

    nt=nt+1;
    sol(nt)=y;
    time(nt)=time_i;

    if (mod(nt,10000)==0)
        fprintf(' %g ',time(nt))
    end

end
cput_fixed=toc;
ufinal_fixed=sol(nt);
uexact=yexact(time(nt),y0);
numerror_fixed=abs(uexact-ufinal_fixed);
fprintf(' time = %g ufinal = %g uexact = %g error_fixed = %g dt_fixed = %g \n',time(nt),ufinal_fixed,uexact,numerror_fixed,dt_fixed);

figure(5)
%subplot(1,3,1)
loglog(dtmax(1:jmax),numerror(1:jmax),'-bo','LineWidth',2,...
    'MarkerSize',14,'MarkerFaceColor',[0,0,1])
hold on
loglog(dtmin(1:jmax),numerror(1:jmax),'-rs','LineWidth',2,...
    'MarkerSize',14,'MarkerFaceColor',[1,0,0])
loglog(dt_fixed,numerror_fixed,'-gs','LineWidth',2,...
    'MarkerSize',14,'MarkerFaceColor',[0,1,0])
axis([min(dtmin)/2,max(dtmax)*2,min(numerror)/2,max(numerror)*2]);
%    dtav=sqrt(dtmax(1:jmax).*dtmin(1:jmax));
dtav=dtmax(1:jmax)/2;
loglog(dtav(1:jmax),numerror(ceil((jmax+1)/2))*(dtav(1:jmax)/dtav(ceil((jmax+1)/2))).^1,'k--')
loglog(dtav(1:jmax),numerror(ceil((jmax+1)/2))*(dtav(1:jmax)/dtav(ceil((jmax+1)/2))).^4,'k--')
loglog(dtav(1:jmax),numerror(ceil((jmax+1)/2))*(dtav(1:jmax)/dtav(ceil((jmax+1)/2))).^5,'k--')
ylabel('error')
xlabel('dt')
hold off
figure(6)
subplot(1,2,1)
loglog(ntj(1:jmax),numerror(1:jmax),'-ro','LineWidth',2,...
    'MarkerSize',14,'MarkerFaceColor',[0,0,1])
xlabel('nt')
ylabel('error')
hold on
loglog(ntj(jmax),numerror_fixed,'-gs','LineWidth',2,...
    'MarkerSize',14,'MarkerFaceColor',[0,1,0])
legend('adaptive','fixed')
hold off
subplot(1,2,2)
loglog(cput(1:jmax),numerror(1:jmax),'-ro','LineWidth',2,...
    'MarkerSize',14,'MarkerFaceColor',[0,0,1])
hold on
loglog(cput_fixed,numerror_fixed,'-gs','LineWidth',2,...
    'MarkerSize',14,'MarkerFaceColor',[0,1,0])
xlabel('cpu time')
ylabel('error')
if (1==2)
    figure(7)
    semilogx(tolj(1:jmax),numerror(1:jmax)./tolj(1:jmax),'ro','LineWidth',2,...
        'MarkerSize',14,'MarkerFaceColor',[0,0,1])
    xlabel('tolerance')
    ylabel('error/tolerance')
end
legend('adaptive','fixed')
hold off

if (adapt==1)
    figure(4)
    hold off
end
if (adapt==0)
    figure(10)
    hold off
end
figure(11)
hold off
if (adapt==1)
    figure(12)
    hold off
    figure(13)
    hold off
    figure(20)
    hold off
    figure(21)
    hold off
    %figure(14)
    %hold off
end
%figure(15)
%hold off

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tout,yout,dt,success,errorestimate]=steprk4opt(time,y,dt,tol,tmax,richardson,adapt)

if (adapt==1)
    p=5;
    pm=p-1;
    [t1,y1]=rk4(time,y,dt/2);
    [t2,y2]=rk4(t1,y1,dt/2);
    [t3,y3]=rk4(time,y,dt);
    diff=y2-y3;
    errorestimate=diff/(1-2^(-pm));
    if (errorestimate==0)
        fprintf(' \n !!! errorestimate =  %g !!! \n', errorestimate)
    end
    dtopt=((tol/tmax)*(dt^p)/abs(errorestimate))^(1/pm);
    %dtopt=((tol*y/tmax)*(dt^p)/abs(errorestimate))^(1/pm);
    if dtopt > dt
        tout=t3;
        if (richardson==0)
            yout=y2;
        else
            yout=2^pm*y2/(2^pm-1)-y3/(2^pm-1);
        end
        %error=yexact(dt,y)-y3;
        %fprintf(' error = %g estimate = %g \n',error,errorestimate);
        success=1;
        %    dt=min(.9*dtopt,2*dt);
        dt=.8*dtopt;
    else
        fprintf(' not accepting time step time = %g  dt = %g dtopt = %g \n',time,dt,dtopt);
        tout=time;
        yout=y;
        dt=.8*dtopt;
        success=0;
    end
else
    [tout,yout]=rk4(time,y,dt);
    success=1;
    errorestimate=1;
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tout,yout,dt,success,errorestimate]=steprk2opt(time,y,dt,tol,tmax,richardson,adapt)

if (adapt==1)
    p=3;
    pm=p-1;
    [t1,y1]=rk2(time,y,dt/2);
    [t2,y2]=rk2(t1,y1,dt/2);
    [t3,y3]=rk2(time,y,dt);
    diff=y2-y3;
    errorestimate=diff/(1-2^(-pm));
    dtopt=((tol/tmax)*(dt^p)/abs(errorestimate))^(1/pm);
    %dtopt=((tol*y/tmax)*(dt^p)/abs(errorestimate))^(1/pm);
    if dtopt > dt
        tout=t3;
        if (richardson==0)
            yout=y2;
        else
            yout=2^pm*y2/(2^pm-1)-y3/(2^pm-1);
        end
        %error=yexact(dt,y)-y3;
        %fprintf(' error = %g estimate = %g \n',error,errorestimate);
        success=1;
        %    dt=min(.9*dtopt,2*dt);
        dt=.9*dtopt;
    else
        fprintf(' not accepting time step time = %g  dt = %g dtopt = %g \n',time,dt,dtopt);
        tout=time;
        yout=y;
        dt=.9*dtopt;
        success=0;
    end
else
    [tout,yout]=rk2(time,y,dt);
    success=1;
    errorestimate=1;
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tout,yout,dt,success,errorestimate]=stepeuleropt(time,y,dt,tol,tmax,richardson,adapt)

if (adapt==1)
    [t1,y1]=euler(time,y,dt/2);
    [t2,y2]=euler(t1,y1,dt/2);
    [t3,y3]=euler(time,y,dt);
    %diff=sqrt((y3(1)-y2(1))^2+(y3(2)-y2(2))^2);
    diff=y2-y3;
    errorestimate=2*diff;
    dtopt=(tol/tmax)*(dt^2)/abs(errorestimate);
    if dtopt > dt
        tout=t3;
        if (richardson==0)
            yout=y2;
        else
            yout=2*y2-y3;
        end
        %error=yexact(dt,y)-y3;
        success=1;
        dt=0.8*dtopt;
        %if (dt>1e-3)
           % fprintf(' error = %g estimate = %g dt = %g tout = %g \n',error,errorestimate,dt,tout);
        %end
    else
        fprintf(' not accepting time step time = %g  dt = %g dtopt = %g \n',time,dt,dtopt);
        tout=time;
        yout=y;
        [time,dt,dtopt];
        dt=.8*dtopt;
        success=0;
    end
else
    [tout,yout]=euler(time,y,dt);
    success=1;
    errorestimate=1;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [time,y]=rk4(time,y,dt)
k1=dt*F(time,y);
k2=dt*F(time+dt/2,y+k1/2);
k3=dt*F(time+dt/2,y+k2/2);
k4=dt*F(time+dt,y+k3);
y=y+(k1+2*k2+2*k3+k4)/6;
time=time+dt;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [time,y]=rk2(time,y,dt)
k1=dt*F(time,y);
k2=dt*F(time+dt,y+k1);
y=y+(k1+k2)/2;
time=time+dt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [time,y]=euler(time,y,dt)
y=y+dt*F(time,y);
time=time+dt;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Fy] = F(x,y)
global ifunction

%pause(.0005);

if (ifunction==1)
    Fy=-1/(4*y^3);
end
if (ifunction==2)
    Fy=y^2;
end
if (ifunction==3)
    Fy=y;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y=yexact(t,y0)
global ifunction

if (ifunction==1)
    y=sqrt(sqrt(y0.^4-t));
end
if (ifunction==2)
    y=1./(1./y0-t);
end
if (ifunction==3)
    y=y0.*exp(t);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tout,yout,dt,success]=stepeuler(time,y,dt,tol)
global diff;
[t1,y1]=euler(time,y,dt/2);
[t2,y2]=euler(t1,y1,dt/2);
[t3,y3]=euler(time,y,dt);
diff=sqrt((y3-y2)^2);
if diff < tol
    tout=t3;
    yout=2*y2-y3;
    success=1;
    if diff < tol/2
        dt=dt*sqrt(2);
    end
else
    tout=time;
    yout=y;
    dt=dt/sqrt(2);
    success=0;
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function marker=make_markers

colors=['b','r','c','m','g','k','y'];
symbols=['o','+','*','x','s','d','h'];

is=0;
for isymbol=1:length(symbols)
    for icolor=1:length(colors)
        is=is+1;
        marker{is}=strcat('-',colors(icolor),symbols(isymbol));
        marker{100+is}=strcat('--',colors(icolor),symbols(isymbol));
    end
end

end
