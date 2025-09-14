%% Givens
data = readmatrix('Case_b1.txt');
L = 4; %length of the beam in m
E = 69 * 10^9; %Youngs modulus (Pa)
I = 2.474 *10^-6; %Second moment of inertia (m^4)
t = (1/16) *0.0254;%Strut outside washer thickness (M)
d = (3/8) * 0.0254; %Strut diameter (M)
A = (pi/4)*(d^2 - (d-t)^2) %Beam area cross section (M^2)


P = data(:,1) .* 4.44822; %conver to N
F0 = data(:,2).* 4.44822; %convert to N
F1 = data(:,3).* 4.44822; %convert to N
F2 = data(:,4).* 4.44822; %convert to N
F3D = ((data(:,5) - data(1,5)).* 4.44822); %convert to N and 0 data compared to first value;
deflection = data(:,6).* 0.0254; %conver to M

%% Averaging all load cases to see where load is 
F0ave = F0(51);
F1ave = F1(51);
F2ave = F2(51);

Freaction1ave = (F0 + F1);
Freaction2ave = F2ave;



%% Solve for a and b:
a = mean(F2(11:100).*(L./P(11:100)));
b0 = mean(Freaction1ave(11:100).*(L./P(11:100)));
Lcheck = a+b0;

b = 4-a;

midspanDefmodel = @(b) (1./(E.*I)) .* ((P.*L^2.*b)./48 + (P.*(b-a))./48 + (P.*b.*(L^2-b^2))./48); %m
midspanModelStress = @(a) (2.526*10^4).*P.*a; %
modelStress = midspanModelStress(a);
midspanModelInternalForce = modelStress .* A; %N

modeledDef = midspanDefmodel(b);
% modeledDefave = mean([modelStress((10*1+1):(10 + 10*1)); modelStress((100 - 10*1 +1):(110 - 10*1))]);

%% For loop to create averages of model and measured forces for each loading condition (excluding 0 load condition)
for i=1:5 %Only 5 times, cause we have 5 variations of loads!

    modeledDefave(i) = mean([modeledDef((10*i+1):(10 + 10*i)); modeledDef((100 - 10*i +1):(110 - 10*i))]);
    modeledInternalForceave(i) = mean([midspanModelInternalForce((10*i+1):(10 + 10*i)); midspanModelInternalForce((100 - 10*i +1):(110 - 10*i))]);

    measuredInternalForceave(i) = mean([F3D((10*i+1):(10 + 10*i)); F3D((100 - 10*i +1):(110 - 10*i))]);
    measuredDefave(i) = mean([deflection((10*i+1):(10 + 10*i)); deflection((100 - 10*i +1):(110 - 10*i))]);

end

deflectionError = abs(modeledDefave - measuredDefave);
internalForceError = abs(modeledInternalForceave - measuredInternalForceave);

comparisonTable = [1000.*modeledDefave', 1000.*measuredDefave',1000.*deflectionError', modeledInternalForceave', measuredInternalForceave', internalForceError']; %comparison table in mm and N
xlswrite('part2task2table.xlsx',comparisonTable);

