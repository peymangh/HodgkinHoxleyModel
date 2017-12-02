%%Generate a spike and Gating variables of the neuron with Hodgkin-Hoxley model of neuron.
%Just write this command on your command window: answer= hodgkinHoxley();
%Then you have a matrix of answer=[t,v,n,m,h] including Time, Voltage of membrane, n variable, m variable,
% and h variable respectively.
%%implementation by Peyman Ghasemi - pe.ghasemi@ut.ac.ir - University of Tehran - Summer 2015

%%H-H
function answer= hodgkinHoxley()
%answer=[t,v,n,m,h]

N=10000;   %Number of Points
a=0;       %Start Time
b=20;      %End Time

%Initial Conditions
V0=0; % Initial Membrane voltage
an= 0.01*((10-V0)/(exp((10-V0)/10)-1));  %Initial alpha n
Bn= 0.125*exp(-V0/80);                   %initial beta n
am= 0.01*((25-V0)/(exp((25-V0)/10)-1));  %and so on
Bm=4*exp(-V0/18);
ah=0.07*exp(-V0/20);
Bh=1/(exp((30-V0)/10)+1);
m0=am/(am+Bm); % Initial m-value
n0=an/(an+Bn); % Initial n-value
h0=ah/(ah+Bh); % Initial h-value
alpha= [V0,n0,m0,h0];  %to send initial values to runge-cutta



%%Runge-Cutta-
m = size(alpha,1);
if m == 1
   alpha = alpha';
end

h = (b-a)/N;        %the step size
t(1) = a;
w(:,1) = alpha;     %initial conditions

for i = 1:N
   k1 = h*f(t(i), w(:,i));
   k2 = h*f(t(i)+h/2, w(:,i)+0.5*k1);
   k3 = h*f(t(i)+h/2, w(:,i)+0.5*k2); 
   k4 = h*f(t(i)+h, w(:,i)+k3);
   w(:,i+1) = w(:,i) + (k1 + 2*k2 + 2*k3 + k4)/6;
   t(i+1) = a + i*h;
end
%---------------------------------------------------

answer=[t;w];   %Output ----> answer=[t,v,n,m,h]

end
%% H-H Equations
function dy= f(t,y) %This function includes Hodgkin-Hoxley Equations

%Applied Current
I=5*heaviside(t-2)- 5*heaviside(t-2.5) + 40*heaviside(t-10) - 40*heaviside(t-10.5);

%Give name to parameters
v=y(1);
n=y(2);
m=y(3);
h=y(4);

%Constants
Ek= -12; %mV
ENa= 120; %mV
EL= 10.6; %mV

gkbar= 36; %mS/cm2
gNabar= 120; %mS/cm2
gL= 0.3;   %mS/cm2

C=1;  %uF/cm2

%Define Alpha & Betha
an=0.01* ( (10-v)/( exp((10-v)/10)-1));
bn=0.125*exp(-v/80);
am=0.1* (25-v)/(exp((25-v)/10)-1);
bm= 4*exp(-v/18);
ah= 0.07*exp(-v/20);
bh=1/(exp((30-v)/10)+1);

%H-H Equations
vdot= (I - gkbar*n^4*(v-Ek) - gNabar*m^3*h*(v-ENa) - gL*(v-EL))/C;
ndot= an*(1-n) - bn*n;
mdot= am*(1-m) - bm*m;
hdot= ah*(1-h) - bh*h;

%Send diffrential of parameters to runge-cutta
dy = [vdot;ndot;mdot;hdot];

end
%%
%Peyman Ghasemi
%pe.ghasemi@ut.ac.ir
%University of Tehran - Iran
%Summer 2015
