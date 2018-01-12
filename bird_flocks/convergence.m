function convergence

clear
format long
y0=1;

dt0=0.1;
tmax=.99;

% forward Euler: method=1
% RK4: method=2
method=2;

choice=3;
switch choice
    % gets to a slope of 2
    case 1
        jmax=20;
        dtfactor=1.2;
        % gets just to the onset of 4thorder
    case 2
        jmax=35;
        dtfactor=1.2;
        % gets just to the noise floor
    case 3
        jmax=43;
        dtfactor=1.2;
    case 4
        jmax=14;
        dtfactor=2;
        
end



% allocate memory to the arrays to speed up the code
ntmax=1000000;
% sol stores position of particle for all times
sol=zeros(ntmax,1);
% for variable time step need to keep track of deltat and t separately
deltat=zeros(ntmax,1);
t=zeros(ntmax,1);
ntj=zeros(ntmax,1);


for j=1:jmax
  dt=dt0/(dtfactor^(j-1));
  dtj(j)=dt;
end

for j=1:jmax
 dt=dtj(j);
 numsteps=floor(tmax/dt)+1;
 x(1:numsteps)=0;
 y(1:size(y0),1:numsteps)=0; 

 y(1,1)=y0;

 cput1=cputime;

 for k=1:numsteps

  if x(k)+dt > tmax
   dt=tmax-x(k);
  end

 if (method==1)
   [Fy] = F(x(k),y(:,k));
   y(:,k+1) = y(:,k)+dt*Fy';
   x(k+1)=x(k)+dt;
 end
 
 if (method==2)
   [Fy] = F(x(k),y(:,k));
   yk1 = dt*Fy';
 
  [Fy] = F(x(k)+dt/2,y(:,k)+yk1/2);
   yk2 = dt*Fy'; 

  [Fy] = F(x(k)+dt/2,y(:,k)+yk2/2);
  yk3 = dt*Fy';

  [Fy] = F(x(k)+dt,y(:,k)+yk3);
  yk4 = dt*Fy';
  
  y(:,k+1)=y(:,k)+(yk1+2*yk2+2*yk3+yk4)/6;
  y(:,k+1)=y(:,k)+(yk1+2*yk2+2*yk3+yk4)/6;
  x(k+1)=x(k)+dt;
 end

 end

 elapsed_time=cputime-cput1;


% for dy/dt=-1/2y
 yexact=sqrt(1-x);
% for dy/dt=-1/4y^3
% yexact=sqrt(sqrt(1-x));
 solution=y(:,numsteps+1);
 solex=yexact(numsteps+1);
 error(j)=abs(y(1,numsteps)-yexact(numsteps));
 fprintf(' solution = %g solex = %g error = %g numsteps = %g \n',solution,solex,error(j),numsteps);
% demonstrate using incorrect timing for exact solution
% error(j)=abs(y(1,numsteps)-yexact(numsteps-1))
 ntj(j)=numsteps;

 figure(1)
%  subplot(1,3,1)
  plot(x(:),y(1,:),'b+');
  xlabel('Time')
  ylabel('y')
  hold on
  plot(x(:),yexact(:),'r*');
  hold off

 figure(2)
 %subplot(1,3,2)
  loglog(dtj(1:j),error(1:j),'-ro','LineWidth',2,...
         'MarkerSize',14,'MarkerFaceColor',[1,0,0]);
  xlabel('Time Step dt')
  ylabel('Error')
  
figure(3)
%subplot(1,3,3)
    loglog(ntj(1:j),error(1:j),'-ro','LineWidth',2,...
         'MarkerSize',14,'MarkerFaceColor',[0,0,1])
  xlabel('nt')
  ylabel('error')

end

end

    %%%%%%%%%%%%%%%%%%
    
function [Fy] = F(x,y)
 Fy(1)=-1/(2*y(1));
end


