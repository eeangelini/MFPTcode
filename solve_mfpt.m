function Y1 = solve_mfpt(k)
%solves double absorbing boundary BVP for mean first passage time
%on-center, symmertric envelope case (i.e. symmetric well in general)
%takes in values of a from 16 to 30, b is fixed at 15


%Erin Angelini, 10.5.17

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load values of k = k_1*W_max for a from 16 to 30, b fixed at 15
%cd 'mfpt mat'
%load('kValsCtrlAR1p06to2.mat')
%cd '../'

%now get corresponding k for input a value and b = 15
%k = kvec(a-15);

X = [linspace(0,pi,10)]; %range of starting angles alpha
solinit = bvpinit(X,[0 0]);

sol = bvp4c(@mpft,@twobc,solinit); %call bvp solver

x = [linspace(0,pi)];%these are alpha values for plotting
y = deval(sol,x);%elavuate the BVP solution numerically

splitProbs = zeros(1,length(x)); %calculate the splitting probability
for i = 1:length(x)
    splitProbs(i) = split(x(i));
end

Y1 = y(1,:); %(mean time)*split
Y2 = Y1./splitProbs; %mean exit time tau
j = length(Y2)/2; %this is approximately the pi/2 point
midval = Y2(j); %tau(pi/2)

%functions for BVP solver
    function dydx=mpft(x,y)
    dydx = zeros(2,1);
    
    dydx(1) = y(2);
    dydx(2) = -split(x)-k*sin(2*x)*y(2);
         
    end

   function res = twobc(ya,yb)
   res = [ ya(1); yb(1)];
   end
   
%splitting probability function  
    function z=split(x)
        N=1000;
        
        s = [x:(pi-x)/N:pi];
        
        f = exp(0.5*k*(1-cos(2*s)));
        
        s2 = [0:pi/N:pi];
        f2 = exp(0.5*k*(1-cos(2*s2)));
        C = trapz(s2,f2);
        
        z = (1/C).*trapz(s,f);
        
    end

end