clear;

% variable declaration
N = 87500; % no of grid points 87500
x = zeros(N,1); % \rho in the manuscript
s = zeros(N,1); % arc length
z = zeros(N,1);
d_eff = zeros(N,1);
theta = zeros(N,1);
sigma = 1; % gaussian standard deviation
Lambda = 0.157; % sine period
a=1;
delta = 0.00016;

% initial values
z(1) = -7; z(2) = -7;
s(1) = -7; s(2) = -7;


for i=3:N
     if ( abs(theta(i-1))< pi/4 ||  abs(theta(i-1))> 3*pi/4)
        z(i) = z(i-1)+delta;
        s(i) = s(i-1) + sqrt((x(i-1)-x(i-2))^2+(z(i)-z(i-1))^2);
	    d_eff(i) = exp(-s(i)^2/(sigma^2)/2)*a*abs(sin(2*pi*s(i)/Lambda));
        theta(i) = 0.5*asin(d_eff(i));
        dx = tan(theta(i));
        x(i) = x(i-1) + delta*dx;
    else
        x(i) = x(i-1)+delta;
        s(i) = s(i-1) + sqrt((x(i)-x(i-1))^2+(z(i-1)-z(i-2))^2);
        d_eff(i) = exp(-s(i)^2/(sigma^2)/2)*a*abs(sin(2*pi*s(i)/Lambda));
        theta(i) = 0.5*asin(d_eff(i));
        dz = tan(theta(i));
        z(i) = z(i-1) +  delta/dz;
     end

end

% PLOT
figure(1);
plot(z(19000:68500),x(19000:68500),'LineWidth',1);
xlabel('x (mm)')
ylabel('z (mm)')

figure(2);
plot(s(19000:67500),d_eff(19000:67500));
xlabel('s (mm)')
ylabel('Relative Nonlinearity')
hold on;