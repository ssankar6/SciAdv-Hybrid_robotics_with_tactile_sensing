function [v,u]=genspikes(A,scalefactor,a,b,c,d)
% % compute an action potential model based on Izhikevich

% inputs
% % A is a data matrix where the rows are the sensors and the columns are the values of the sensors at different points in time
% % scalefactor is a gain to scale sensor output into excitation currents
% % a,b,c,d are model parameters

% output 
% % v is the membrane potential for the data array
% % u is the recovery variable
spkv = [];
v=-65*ones(size(A,1),size(A,2)+1);    % Initial values of v
u=b.*v;                 % Initial values of u
for t=2:size(A,2)+1           % shift values by 1 since computation depends on initial condition
    I=A(:,t-1)*scalefactor;
    v(:,t)=v(:,t-1)+(0.04*v(:,t-1).^2+5*v(:,t-1)+140-u(:,t)+I);
    u(:,t)=u(:,t-1)+a.*(b.*v(:,t)-u(:,t-1));
    for i=1:size(A,1) % scan each sensor for a spike
        if v(i,t)>=30
          v(i,t-1)=30;
          v(i,t)=c;
          u(i,t)=u(i,t)+d;
        end
    end
end
end