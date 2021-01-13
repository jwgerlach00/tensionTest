% Jacob Gerlach
% jwgerlac@ncsu.edu
% 10/9/2020
% stressStrain_GERLACH.m
%
% Calculates the stress/strain curve for a sample given tension test data.

clear
clc
close all

%% Declarations
sample = xlsread('brass.xlsx'); % excel data
displacement = sample(:, 1); % sample displacement (m)
force = sample(:, 2); % tensile force (N)
length = 38e-3; % sample gage length (m)
radius = 1.5e-3; % sample cross-sectional radius (m)
yOffset = 0.002; % offset for yield
smoothSpan = 25;
strainGuess = 0.100000000000000; % strain to interpolate stress

%% Anonymous Functions
sigma = @(f, r) f/(pi*(r^2)); % stress
epsilon = @(delL, L) delL/L; % strain
MPa = @(Pa) Pa/(10^6); % Pa to MPa conv.
GPa = @(Pa) Pa/(10^9); % Pa to GPa conv.

%% Calculations
force = smooth(force, smoothSpan); % smooth graph
stress = sigma(force, radius);
strain = epsilon(displacement, length);

for k = 1:(size(stress) - 2)
    r(k) = corr(strain(1:(k + 2)), stress(1:(k + 2))); % correlation coeff
end

[~, prop] = max(r); % proportional index

stress_p = stress(prop); % proportional stress (Pa)
strain_p = strain(prop); % proportional strain

eMod = polyfit(strain(1:prop), stress(1:prop), 1); % elastic modulus line
eMod = eMod(1); % elastic modulus slope

yieldLine = eMod*(strain - yOffset); % intersect for yield

[~, yield] = min(abs(stress - yieldLine)); % yield point
stress_y = stress(yield); % yield stress (Pa)
strain_y = strain(yield);

[stress_u, ultimate] = max(stress); % ultimate stress (Pa)
strain_u = strain(ultimate);

stress_f = stress(end); % fracture stress (Pa)
strain_f = strain(end);

stressGuess = interp1(strain, stress, strainGuess); % interpolate stress

resiliance = trapz(strain(1:yield), MPa(stress(1:yield))); % (MJ/m^3)
toughness = trapz(strain, MPa(stress)); % (MJ/m^3)

%% Output
plot(strain, MPa(stress)); % plot
xlabel('Strain');
ylabel('Stress (MPa)');
title('Tension Test');
hold on
toughnessArea = area(strain, MPa(stress)); % area
toughnessArea.FaceColor = 'c';
resilianceArea = area(strain(1:yield), MPa(stress(1:yield)));
resilianceArea.FaceColor = 'r';
plot(strain_p, MPa(stress_p), 'mo', strain_y, MPa(stress_y), 'mo',...
    strain_u, MPa(stress_u), ' mo', strain_f, MPa(stress_f), 'mo');
text(strain_p, MPa(stress_p), '  Proportional');
text(strain_y, MPa(stress_y), '  Yield');
text(strain_u, MPa(stress_u), '  Ultimate');
text(strain_f, MPa(stress_f), '  Fracture');
hold off

fprintf('Proportional stress: %.2f MPa, strain: %.4f\n', MPa(stress_p),...
    strain_p)
fprintf('Yield stress: %.2f MPa, strain: %.4f\n', MPa(stress_y), strain_y)
fprintf('Ultimate stress: %.2f MPa, strain: %.4f\n', MPa(stress_u),...
    strain_u)
fprintf('Fracture stress: %.2f MPa, strain: %.4f\n', MPa(stress_f),...
    strain_f)
fprintf('Elastic modulus: %.2f GPa\n', GPa(eMod));
fprintf('Interpolated stress: %.2f MPa, strain: %.4f\n',...
    MPa(stressGuess), strainGuess);
fprintf('Resiliance: %.2f MJ/m^3\n', resiliance);
fprintf('Toughness: %.2f MJ/m^3\n', toughness);