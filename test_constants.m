global tauHE; % Duration (in days) of host egg stage
global tauHL; % Duration (in days) of host larval stage
global tauHP; % Duration (in days) of host pupal stage
global tauHA; % Duration (in days) of host adualt stage
global tauPL; % Duration (in days) of parasitoid egg and larval stage
global tauPA; % Duration (in days) of parasitoid adult stage
global bHA; % Daily host adult fecundity (number of eggs)
global c; % Competition mortality coefficient
global k; % Interference parameter
global a; % Parasitoid attack rate
global eta; % Encapsulation probability
global dHE; % Background host egg mortality
global dHL; % Background host larval mortality
global dHP; % Background host pupal mortality
global dHA; % Background host adualt mortality
global dPL; % Background parasitoid juvenile mortality
global dPA; % Background parasitoid adult mortality
global eps; % Proportion of mutants
global cH; % Cost of resistance
global cP; % cost of virulence
global A; % Slope of encapsulation trade-off curve
global x; % Background resistance
global y; % Background virulence
global b0; % Maximum host fecunidty
global dPL0; % Minimum parasitoid larval mortality
global N; % Number of host strains
global M; % Number of parasitoid strains

tauHE = 4.3; % Duration (in days) of host egg stage
tauHL = 25.0; % Duration (in days) of host larval stage
tauHP = 7.0; % Duration (in days) of host pupal stage
tauHA = 5.5; % Duration (in days) of host adualt stage
tauPL = 20.0; % Duration (in days) of parasitoid egg and larval stage
tauPA = 2.0; % Duration (in days) of parasitoid adult stage
bHA = 9.4; % Daily host adult fecundity (number of eggs)
c = 5*10^-5; % Competition mortality coefficient
k = 0.01; % Interference parameter
a = 0.01; % Parasitoid attack rate
dHE = 0.0; % Background host egg mortality
dHL = 0.0; % Background host larval mortality
dHP = 0.0; % Background host pupal mortality
dHA = 1.0; % Background host adualt mortality
dPL = 0.1; % Background parasitoid juvenile mortality
dPA = 0.1; % Background parasitoid adult mortality
eps = 0.00001; % Proportion of mutants
A = 0.5; % Slope of encapsulation trade-off curve
x = 10.0; % Background resistance
y = 10.0; % Background virulence
b0 = 80.0; % Maximum host fecunidty
dPL0 = 0.01; % Minimum parasitoid larval mortality
N = 1; % Number of host strains
M = 1; % Number of parasitoid strains
