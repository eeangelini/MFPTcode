function solve_SAB
%numerically evaluates the integral given by the single absorbing boundary
%problem for mean first passage time

%Erin Angelini, 5.16.18

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a=0.001; b=pi/2;
alpha = linspace(a,b);

T = zeros(1,length(alpha)); %vector for mean time solution T

fxn = @(y) (b-sin(b)*cos(b)+sin(y).*cos(y)-y)./(1-cos(2.*y));

for i = 1:length(alpha)
    ymin = a; ymax = alpha(i);
    T(i) = integral(fxn,ymin,ymax);
end

T = 1.*T;

plot(alpha,T,'k-','Linewidth',4)
xlim([0 max(alpha)])
xlabel('Starting angle \alpha');
ylabel('mean time \tau_0(\alpha)');
end