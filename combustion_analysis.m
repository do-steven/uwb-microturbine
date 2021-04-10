%fuel information (hydrocarbon only)
a = 1; %number of carbon atoms
b = 4; %number of hydrogen atoms
xs = 0; %percent-excess
xs_f = 1 + xs/100; %multiplication factor for excess

%initialize variables
syms x1...%moles of CH4
     x2...%moles of O2
     x3...%moles of N2
     x4...%moles of CO2
     x5...%moles of H20
     x6   %moles of N2

%chemical reaction without excess
carbon = a*x1 - x4 == 0; %number of carbon atoms
hydrogen = b*x1 - 2*x5 == 0; %number of hydrogen atoms
oxygen = 2*x2 - 2*x4 - x5 == 0; %number of oxygen atoms
nitrogen = 2*x3 - 2*x6 == 0; %number of nitrogen atoms
air = x3 == (79/21)*x2; %oxygen-nitrogen ratio
fuel = x1 == 1; %per mole of fuel

%initialize system of equations
[A,B] = equationsToMatrix([carbon,hydrogen,oxygen,nitrogen,air,fuel],...
    [x1,x2,x3,x4,x5,x6]);

%balancing chemical equation for STOICHIOMETRIC air
x = [x1;x2;x3;x4;x5;x6];
moles = round(linsolve(A,B),5,'significant'); %balanced chemical equation

%chemical equation for stoichiometric air
moles_stoich = moles;

%air-fuel ratio for stoichiometric air
AF_stoich = round((moles_stoich(2)*32 + moles_stoich(3)*28.01) / ...
    (moles_stoich(1)*(12.01*a + 1.01*b)),5,'significant'); 

%balancing chemical equation if there is EXCESS air
if xs ~= 0
    syms x7 %moles of unreacted O2
    %relationship for excess air and oxygen atoms
    oxygen = xs_f*2*moles(2) - 2*moles(4) - moles(5) - 2*x7 == 0;
    x = cat(1,x,x7);
    x7 = solve(oxygen,x7); %calculate excess O2
    moles(2:3) = moles(2:3)*xs_f; %adjust reactants for excess air
    moles(6) = moles(6)*xs_f; %adjust product for excess air
    chem_eqn_excess = cat(2,x,moles); %balanced chemical equation
end

%chemical equation for excess air
moles_actual = moles;

%air fuel ratio for excess air
AF_actual = round((moles_actual(2)*32 + moles_actual(3)*28.01) /...
    (moles_actual(1)*(12.01*a + 1.01*b)),5,'significant');

%equivalence ratio
phi = round(AF_stoich / AF_actual,5,'significant');

%inlet conditions
T1 = 300;

%initialize Cantera
%initialize reactant mixture
reactant = GRI30;
nsp = nSpecies(reactant);

%defining reactants
ch4 = speciesIndex(reactant,'CH4');%fuel
o2 = speciesIndex(reactant,'O2');%air
n2 = speciesIndex(reactant,'N2');%air

x = zeros(nsp,1);
x(ch4,1) = moles(1);%moles of CH4
x(o2,1) = moles(2);%moles of O2
x(n2,1) = moles(3);%moles of N2
set(reactant,'T',T1,'MoleFractions',x);%defining state of reactants
product = equilibrate(reactant,'HP');%solves for chemical equilibrium

T_ad = round(temperature(product),0,'decimals'); %adiabatic flame temp

