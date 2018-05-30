function BVPplot(Y)
%this function takes a 6x100 vector Y (either tau or splitting probability
%from the BVP solvers solve_mfpt or solve_mfpt_new) and plots it

%Erin Angelini, 4.23.18

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = linspace(0,pi);
plot(x,Y,'-','Linewidth',4)
xlim([0 max(x)])
xticks([0 pi/4 pi/2 3*pi/4 pi])
xticklabels({'0','\pi/4','\pi/2','3\pi/4', '\pi'})
set(gca,'FontSize',30)
xlabel('Starting angle \alpha');
%ylabel('mean time \tau_0(\alpha)');
ylabel('solution T(\alpha) = \pi_0(\alpha)\tau_0(\alpha)')
%ylabel('splitting probability \pi_0(\alpha)');
end
