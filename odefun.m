% Initial Coefficient Values %
C_p = 19.21712; % Recruitment Potential
N = 43250000; # Total Population of Voters
tau_i = 0.00945; Average time spent as a R.A (in years)
tau_a = 10; % Average time spend as a non-R.A (in years)
alpha = 0.15; % Leaving rate
g = 0.2; % Fraction recruited
%f = (y(2)+y(3))/(y(2)+y(3)+y(4)) % fraction of recruiting and non-recruiting activists

% Initial Population Values %
S0 = 11294000; %initial susceptibles 
I0 = 101; %inital recruiters
A0 = 49899; %initial non-recruiting activists
M0 = 216000; %intial inactive party members

% Defining our system of equations %
f=@(t,x)[-0.000047019*x(1)*x(2)+0.15*x(4); % xa(:,1) = S
         0.000009404*x(1)*x(2)-69.44*x(2); % xa(:,2) = I = Number of Recruiting Activists
         0.000037615*x(1)*x(2)*((x(2)+x(3))/(x(2)+x(3)+x(4)))+69.44*x(2)-0.1*x(3); % A = xa(:,3)
         0.000037615*x(1)*x(2)*((x(4))/(x(2)+x(3)+x(4)))+0.1*x(3)-0.15*x(4); % M = xa(:,4)
         ]; 
[t,xa] = ode45(f,[1994, 2000],[S0 I0 A0 M0]); % inputting initial values
   
 % Plotting Results %
% Population of all member types for the Party %
figure(1);
plot(t,(xa(:,2)+xa(:,3)+xa(:,4))/1000000,t,(xa(:,4))/1000000,t,(xa(:,3))/1000000,t,(xa(:,2)+xa(:,3))/1000000) % graph with all values (Party Members, Inactive Members, Activists, Recruiting Activists)
legend({'Total Party Members','Inactive Members', 'Non-Recruiting Activists', 'Activists'})
title('Party Populations');
xlabel('Year');
ylabel('Population (millions)');

% Thresholds Graph %
figure(2);
horzline = 3.8434*ones(100,1);
slopedline = transpose(linspace(3.8434,3.7777,137));
myline=vertcat(line,slopedline);
plot(t,myline,t,43250000./xa(:,1),'--or',t,3.7413,'--*')
ylim([3.7,4]);

%fplot(@(x) 3.86,[1994,1998])
%hold on
%fplot(@(x) (4-0.5*x)/t,[1998,2004]) 
%hold off;

legend({'Rp','Repi','Rext'})
title('Thresholds');
xlabel('Year');
ylabel('persons/person');

%dydt = zeros(4,1);
% S(susceptibles)
%dydt(1) = -(C_p/tau_i*N)*y(1)*y(2) + alpha*y(4);
% I (infectious/recruiting activists)
%dydt(2) = g*(C_p/tau_i*N)*y(1)*y(2) - y(2)/tau_i;
% A (non-recruiting activists)
%dydt(3)= (y(2)+y(3))/(y(2)+y(3)+y(4))*(1-g)*(C_p/tau_i*N)*y(1)*y(2) + y(2)/tau_i - y(3)/tau_a;
% M (inactive members)
%dydt(4)= (1-(y(2)+y(3))/(y(2)+y(3)+y(4)))*(1-g)*(C_p/tau_i*N)*y(1)*y(2) + y(3)/tau_a - alpha*y(4);

%[t,y] = ode45(@(t,y) odefcn(t,y,C_p,N,tau_i,tau_a,alpha,g), tspan, y0);

