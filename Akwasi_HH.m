%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Template for Hodgkin Huxley model of action potential
% in a giant squid axon
% Complete the code by inserting appropriate expressions 
% in the places marked 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all

% Define constants here
Vrest = -60 ; % resting potential, mV
ENa = 54.2 ; % sodium Nernst potential, mV corrected from 52.4
gNa_max = 120  ; % sodium saturation conductance, mS/cm2
%gNa_max = 120*0.3
%% 1) Define other constant related to IK and IL here
%%
EK = -74.7 %potassium Nernst potential, mV
gK_max = 36 %potassium saturation conductance, mS/cm2
gK_max = 36 *0.5
EL = -43.256 %leak Nernst potential, calculated using the goldmann
gL_max = 0.3 % leakage saturation conductance, mS/cm2
Cm = 1 %membrane capacitance 1uF/cm2

% time stepping and stimulus related constants
deltaT = 0.001 ; % time step, millisec
tStart = -1.000 ; % start time, millisec
tEnd = 15.000 ; % end time, millisec
nStep = ceil((tEnd-tStart)/deltaT) % number of time steps
outputInterval = 20 ; % number of time steps between screen output
StimDur=2; % duration of the stimulus, millisec
Jstim=150; % strength of stimulus


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set initial value of state variables
Vm = Vrest ; % membrane potential, mV
m = 0.05293 ; % initial value of m gate
%% 2) similarly, declare initial values of h and n gates here
%%

h = 0.59612; %initial value for h gate
n = 0.31768; %initial value of n gate


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Storage for variables to plot
plot_Vm = zeros(nStep,1) ;
plot_m = zeros(nStep, 1);
plot_h = zeros(nStep, 1);
plot_n = zeros(nStep, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print a heading for the screen output
display('Hodgkin-Huxley squid axon model');
fprintf('\t Time \t Vm \t n \t\t m \t\t h\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start the simulation
tNow = tStart ;
for iStep = 1:nStep
% Compute ion currents & stimulus current
JNa = gNa_max*m*m*m*h*(Vm-ENa) ;
%% 3) write similar expressions for JK and IL here
%%
%JK = gK_max*n*n*n*n*(Vm - EK);
JK = gK_max*n*n*n*n*(Vm - EK);
JL = gL_max*(Vm-EL);

if( 0<=tNow && tNow<StimDur ) % apply stimulus current at t=0
Jm = Jstim ;
else
Jm = 0 ;
end

% Compute gates' opening and closing rates

alpha_m = 0.1*(25-Vm)/(exp(0.1*(25-Vm))-1);
beta_m = 4*exp(-Vm/18);
m = m + deltaT*(alpha_m*(1-m)-beta_m*m);

%% 4) Write similar expressions for n and h gates here
%%
alpha_n = (0.01*(10-Vm))/(exp((10-Vm)/10)-1);
beta_n = 0.125*exp(-Vm/80);
n = n + deltaT*(alpha_n*(1-n) - beta_n*n);

alpha_h = 0.07*exp(-Vm/20);
beta_h = 1/(exp((30-Vm)/10)+1);
h = h + deltaT*(alpha_h*(1-h) - beta_h*n);

% Compute change in state variables
deltaVm = -deltaT * (JNa + JK + JL - Jm)/Cm ;

% Record/display state variables & other values
plot_Vm(iStep) = Vm ;
plot_m(iStep) = m;
plot_n(iStep) = n;
plot_h(iStep) = h;
plot_time(iStep) = tNow;
if mod(iStep,outputInterval) == 0
fprintf('%8.2f %7.3f %7.5f %7.5f %7.5f\n', ...
tNow, Vm, n, m, h) ;
end % if mod(tNow)

% Update state variables
Vm = Vm + deltaVm ; % new Vm = current Vm + change in Vm
tNow = tStart + iStep*deltaT ; % increment the time
end % for iStep

DataNa30(:,1) = plot_time
DataNa30(:,2) = plot_Vm
% Plot the gates, probabilities, currents, Vm, etc
% subplot(2,1,1)
plot(plot_time, plot_Vm); grid on ;
xlabel('time (ms)')
ylabel('Transmembrane Voltage(Vm), mV')
title('Time course of Action Potential')

%% 5) Similarly plot gates n, m and h against time 
%%
figure
plot(plot_time, plot_m); grid on ;
hold on
plot(plot_time, plot_h);
hold on
plot(plot_time, plot_n);
legend('m', 'h', 'n')
xlabel('time (ms)')
ylabel('gating activation/ inactivation probability')
title('gating probabilities (m, h, n) vrs time')

% end of file