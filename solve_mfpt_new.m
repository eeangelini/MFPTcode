function splitProbs = solve_mfpt_new(k,j)
%solves double absorbing boundary BVP for mean first passage time
%for off center, asymmetric envelope (i.e. asymmetric well in general)
%k is Wmax, j is Wmin, obtain these from Main.m file for each AR

%Erin Angelini, 3.6.18

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


solinit = bvpinit(linspace(0,pi,10),[0 0]);

sol = bvp4c(@mpft,@twobc,solinit);

x = linspace(0,pi);
y = deval(sol,x);

splitProbs = zeros(1,length(x));
for i = 1:length(x)
    splitProbs(i) = split(x(i));
end

Y1 = y(1,:); %tau*split
Y2 = Y1./splitProbs; %mean exit time tau
%J = length(Y2)/2;
%midval = Y2(J);


    function dydx=mpft(x,y)
    dydx = zeros(2,1);
    
    dydx(1) = y(2);
    dydx(2) = -split(x)-0.001*(k*cos(x)-j*sin(x)*0.5*(1+tanh(20*(x-pi/2)))+j*cos(x)*0.5*(1-tanh(20*(x-pi/2))^2))*y(2);
         
    end

   function res = twobc(ya,yb)
   res = [ ya(1); yb(1)];
   end 
   
   
    function z=split(x)
        N=1000;
        
        s = [x:(pi-x)/N:pi];
        
        f = exp(0.001*k*sin(s)-0.001.*j.*cos(s).*0.5.*(1+tanh(20.*(s-pi/2))));
        
        s2 = [0:pi/N:pi];
        f2 = exp(0.001.*k.*sin(s2)-0.001.*j.*cos(s2).*0.5.*(1+tanh(20.*(s2-pi/2))));
        C = trapz(s2,f2);
        
        z = (1/C).*trapz(s,f);
        
    end

end