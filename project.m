%% demo ODE solution. harmonic oscillator
close all; clear;
%% Define the ODE parameters
k2m = 2;
B = 0.1;
w = sqrt(k2m);
T = 2*pi./w;
% x' = v
% v' = -k*x - B*v
xp = @(t, v) v;
vp = @(t, x, v) -k2m*x - B*v;
ta = 0;
tb = 100;
% x = A0*sin(w*t+phi0); -> A0*sin(phi0) = x0
% v = A0*w*cos(w*t+phi0); -> A0*w*cos(phi0) = v0
% phi0 = arctan(x0/v0*w)
% A0 = x0/sin(phi0);
x0 = 0; %
v0 = 1;
% analytical solution:
phi0 = atan(x0./v0.*w);
A0 = sqrt(x0.^2 + (v0./w).^2);
disp(['Analytical solution:x=',num2str(A0,'%.2f'),'*sin(',num2str(w,'%.2f'),'*t+',num2str(phi0,'%.2f'),')']);
for npts = 1001
    tt = linspace(ta,tb,npts);
    h = (tb-ta)/(npts-1);
    % theoretical solution
    x_th = A0.*sin(w.*tt+phi0);
    v_th = w.*A0.*cos(w.*tt+phi0);
    % % modified Euler's method.
    % initialize the solution grid
    x_me = zeros(1,npts);
    v_me = zeros(1,npts);
    % Set the initial/boundary condition.
    x_me(1) = x0;
    v_me(1) = v0;
    for i = 1:npts-1
        % first get y' at half step
        % xp_ct = x_me(i) + h/2*xp(tt(i), v_me(i));
        % vp_ct = v_me(i) + h/2*vp(tt(i), x_me(i)); 
       
        % xp_ct = v_me(i) - k2m.*x_me(i).*(h/2); 
        % vp_ct = -k2m.*(x-me(i)+v_me(i).*(h/2)) - ;
        x_ct = x_me(i) + (h./2)*xp(tt(i), v_me(i));
        v_ct = v_me(i) + (h./2)*vp(tt(i), x_me(i), v_me(i));
        % Euler's update
        x_me(i+1) = x_me(i) + h*xp(tt(i)+h/2, v_ct);
        v_me(i+1) = v_me(i) + h*vp(tt(i)+h/2, x_ct, v_ct);
    end
end
v_fwd = zeros(1,length(x_me));
v_bwd = zeros(1,length(x_me));
v_ct = zeros(1,length(x_me));
v_bwd(1) = v0;
v_ct(1) = v0;
v_fwd(length(v_fwd)) = 0;
v_ct(length(v_fwd)) = 0;
for i = 1:length(x_me)
    if(i<length(x_me))
        disp(i + "forwards")
        v_fwd(i) = (x_me(i+1) - x_me(i))/h;
    end
    if(i>1)
        disp(i + "backwards")
        v_bwd(i) = (x_me(i) - x_me(i-1))/h;
    end
    if(i>1 && i<length(x_me))
        disp(i + "centered")
        v_ct(i) = (x_me(i+1) - x_me(i-1))/(2*h);
    end

end


%% Plot the results of the Euler's method
lwid = 2; font = 16;
nrow = 2; ncol = 1;
figure;
subplot(nrow,ncol,1);
hold on; box on; grid on; axis tight;
set(gca,'LineWidth',lwid,'Fontsize',font);
% plot(tt, x_th, '-','LineWidth',lwid, 'DisplayName','Theoretical');
plot(tt, x_me, '-.','LineWidth',lwid, 'DisplayName','Modified Euler');
xlabel('t');
ylabel('x');
title(['Solution of harmonic oscillator using modified Euler''s Method with',num2str(npts-1,'%d'),' slices']);
subplot(nrow,ncol,2);
hold on; box on; grid on; axis tight;
set(gca,'LineWidth',lwid,'Fontsize',font);
% plot(tt, v_th, '-','LineWidth',lwid, 'DisplayName','Theoretical');
plot(tt, v_me, '-.','LineWidth',lwid, 'DisplayName','Modified Euler');
xlabel('t');
ylabel('v');
hold on; box on; grid on; axis tight;
legend show;
hold on; box on; grid on; axis tight;
set(gca,'LineWidth',lwid,'Fontsize',font);
figure();
% plot(tt, v_th, '-','LineWidth',lwid, 'DisplayName','Velocity Finite Difference Method Approximation');
plot(tt, v_fwd, '-.','LineWidth',lwid, 'DisplayName','Forward Finite Difference Method');
hold on;
plot(tt, v_bwd, '-.','LineWidth',lwid, 'DisplayName','Backward Finite Difference Method');
hold on;
plot(tt, v_ct, '-.','LineWidth',lwid, 'DisplayName','Centered Finite Difference Method');
xlabel('t');
ylabel('v');
hold on; box on; grid on; axis tight;
legend show;
