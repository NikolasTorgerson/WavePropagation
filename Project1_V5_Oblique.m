%% Wave Propigation Through Multilayers
%% Nikolas Torgerson, Gordon Piegan
%% EELE 432

close all; clear all; clc;

%Read Me
%The following code is divided into several parts, just below this are
%three examples, the first being a Distributed Bragg Reflector (also known
%as a dielectric mirror), the second is a problem pulled directly from
%example 5.9 within the book. The third example is meant to be a framework
%that you can add your own specifications into.
%Following the examples is a section that uses the entered information to
%calculate the impedence of each layer, as well as the propigation, angles
%and reflection and transmission coefficients of each interface. This
%information is then used in the final matrix form and the reflected and
%transmitted components of the entire structure are calculated. In
%addition, this code iterates through this process for a specificed range
%of wavelengths, this means that you must calculate the wavelength of your
%desired frequency or at least know what section of the data you would like
%to look at.
%When entering the epsilon (permitivity) of the material, the first and
%last layer represent the medium the wave starts in and ends in and are assumed
%to be infinite. Likewise,the first and last layer thickenss do not actually
%matter as they are assumed to be infinite.

%Initialize all values
%--------------------------------------------------------------------------
u0 = 4*pi*10^(-7);
e0 = 8.85*10^(-12);

syms Eto_II Ero_II
syms Eto_T Ero_T
syms Eio eqn freq

Eio = 1;

Example_to_execute = 1;

switch Example_to_execute
    case 1
% Example 1 Distributed Bragg Reflector-----------------------------------
% This example seeks to simulate a Distributed Bragg Reflector which is a
% structure composed of alternating layer pairs. This structure can also be
% referred to as a dielectric mirror. A design frequency is chosen for the
% structure and the layer thicknesses are calculated to be a quarter of the
% wavelength in the material. The frequency of the wave is varied and the
% resulting reflection and transmission coefficients are graphed.

N = 13; %Number of layers of Bragg Reflector
Eps1 = 8; Eps2 = 4; %Epsilon of the first and second layers respectivly, the first layer must be higher than the second
Mu1 = 1; Mu2 = 1; %Mu values of the layers
freqDesign = 0.8*10^15; %Design frequency for the structure

%Create the array that will store the information about the layers
er_layers = zeros(N*2+2, 0);
ur_layers = zeros(N*2+2, 0);
er_layers(1) = 1;
ur_layers(1) = 1;
er_layers(N*2+2) = 1;
ur_layers(N*2+2) = 1;
for val = 1:N
    er_layers(val*2) = Eps1;
    er_layers(1+val*2) = Eps2;
    ur_layers(val*2) = Mu1;
    ur_layers(1+val*2) = Mu2;
end

wavelengths = zeros(N*2+2, 0); %Calculate the wavelength for each layer
for val = 1:size(wavelengths)
    wavelengths(val) = 1/(freqDesign*sqrt(er_layers(val)*e0*ur_layers(val)*u0));
end

qwavelenghts = wavelengths/4; %Divide the wavelengths by 4
layer_thickness = qwavelenghts; 
incident_angle = 89.9999;

range = [200:1:650]*10^(-9); %This is the range of wavelengths that will be displayed

%End of example 1 --------------------------------------------------------

    case 2
%Example 2-----------------------------------------------------------------
% This is from example 5-9 in the book. The specifications being that
% epsilon of the middle layer is 2.56, the wave is normally incident, and
% you are supposed to find the thickness of the middle layer that will make
% the reflection coefficient go to 0. The value is calculated to be half
% integers of 1.875cm and it can be seen that when the thickness is set to
% 0.0281 (1.5 x 1.875cm), the reflection coefficient goes to 0.

er_layers       = [1 2.56 1];    %Relative epsilon values for each material
ur_layers       = [1 1 1];    %Relative mu values for each material
layer_thickness = [1 1 1];    %Thickness of each material layer (meters)
incident_angle  = 0;           %Angle of incident angle (degrees)
freqDesign      = 10*10^9;      %Frequency of wave (Hz)
Eio             = 1;            %Initial amplitude of incoming wave (V/m)
layer_thickness(2) = (3/2)/(freqDesign*sqrt(er_layers(2)*e0*ur_layers(2)*u0));
wavelength_range = 1/(freqDesign*sqrt(u0*e0))
range = [200:1:650]*wavelength_range*0.005;
%End of Example 2----------------------------------------------------------

    case 3
%Modifiable ---------------------------------------------------------------
er_layers       = [1 1.38 1.701 1.52];    %Relative epsilon values for each material
ur_layers       = [1 1 1 1];    %Relative mu values for each material
layer_thickness = [1 1 1 1];    %Thickness of each material layer (meters)
incident_angle  = 45;           %Angle of incident angle (degrees)
freqDesign       = 10*10^9;      %Frequency of wave (Hz)
Eio             = 1;            %Initial amplitude of incoming wave (V/m)
wavelength_range = 1/(freqDesign*sqrt(u0*e0))
range = [80:1:400]*wavelength_range*0.005;
wavelengths = [1 1 1 1]; %Calculate the wavelength for each layer
for val = 1:length(wavelengths)
    wavelengths(val) = 1/(freqDesign*sqrt(er_layers(val)*e0*ur_layers(val)*u0));
end
layer_thickness = wavelengths/4; %Quarter Wavelengths
%End of Example 3----------------------------------------------------------

    case 4
%Example 4-----------------------------------------------------------------
er_layers       = [1 1.896129 1 2.127 1];    %Relative epsilon values for each material
ur_layers       = [1 1 1 1 1];    %Relative mu values for each material
layer_thickness = [1 45*10^(-9) 301*10^(-9) 56*10^(-9) 1];    %Thickness of each material layer (meters)
incident_angle  = 0;           %Angle of incident angle (degrees)
Eio             = 1;            %Initial amplitude of incoming wave (V/m)
range = [200 : 1 : 1000]*10^(-9);
%End of Example 4----------------------------------------------------------
end

%--------------------------------------------------------------------------
wavelengths_to_test = length(range);


ref_data_reflection_S    = zeros(1, wavelengths_to_test);
ref_data_transmission_S  = zeros(1, wavelengths_to_test);
ref_data_reflection_P   = zeros(1, wavelengths_to_test);
ref_data_transmission_P = zeros(1, wavelengths_to_test);
totals_S = zeros(1, wavelengths_to_test);
totals_P = zeros(1, wavelengths_to_test);

for wavelength = 1:wavelengths_to_test %Beginning of for loop to test for freqencies in range
    frequency = 1/(range(wavelength)*sqrt(u0*e0));
 
%--------------------------------------------------------------------------
%The following declarations prepare several matricies for holding the
%information that will be calculated based on the user's input.
materials_number = length(er_layers);
interface_number = materials_number-1;

impedence_layers    = zeros(1, materials_number);   %Array to hold the calculated impedences of each material
k_layers            = zeros(1, materials_number);   %Array to hold propigation coefficients for each material
k_x                 = 0;                        %Propigation constant in the x direction. This value is the same for each layer
k_z_layers          = zeros(1, materials_number);   %Propigation constants in the z directions for each layer. This value changes.
angles              = zeros(1, materials_number);   %Angle of wave in each layer
angles(1)           = incident_angle*(pi/180);  %First angle is equal to the incident angle

reflect_coeff_P    = zeros(1, interface_number); %Parallel reflection coefficients for each interface
transmit_coeff_P   = zeros(1, interface_number); %Parallel transmission coefficients for each interface
reflect_coeff_S     = zeros(1, interface_number); %Perpendicular reflection coefficients for each interface
transmit_coeff_S    = zeros(1, interface_number); %Perpendicular transmission coefficients for each interface

%--------------------------------------------------------------------------
%Iterative method for calculating the impedences, propigation constants(k),
%and angles for each material layer in the structure.
for m = 1:length(er_layers)
    impedence_layers(m) = sqrt((ur_layers(m)*u0)/...    %Calculates the impedences in each layer
                               (er_layers(m)*e0)); 
    %Calculates the propigation constant for each layer
    k_layers(m) = frequency*2*pi*sqrt(e0*er_layers(m)*u0*ur_layers(m));
    %Calculates the kx propigation coefficient (due to boundary conditions 
    %this is the same for every layer in the structure
    k_x = k_layers(1)*sin(angles(1));
    %Calculates the kz propigation coefficient based off of k and kx
    k_z_layers(m) = sqrt(k_layers(m)^2 - k_x^2);
end
for m = 2:length(angles)
    %Calculates the angle using Snell's law of refraction for each layer
    %based off of the initial incident angle and the propigation constant k
    %of each layer.
    angles(m) = asin(sin(angles(m-1))*k_layers(m-1)/k_layers(m));
end
%--------------------------------------------------------------------------
%Iterative method for calculating the reflection and transmission
%coefficients for the parallel and perpendicular components utilizing the
%equations given by the book.
for m = 1:(length(er_layers)-1)
    reflect_coeff_P(m) = (impedence_layers(m+1)*cos(angles(m+1))-impedence_layers(m)*cos(angles(m)))/...
                          (impedence_layers(m)*cos(angles(m))+impedence_layers(m+1)*cos(angles(m+1)));
    transmit_coeff_P(m)= (2*impedence_layers(m)*cos(angles(m)))/...
                          (impedence_layers(m)*cos(angles(m))+impedence_layers(m+1)*cos(angles(m+1)));
    reflect_coeff_S(m)  = (impedence_layers(m+1)*cos(angles(m))-impedence_layers(m)*cos(angles(m+1)))/...
                          (impedence_layers(m+1)*cos(angles(m))+impedence_layers(m)*cos(angles(m+1)));
    transmit_coeff_S(m) = (2*impedence_layers(m+1)*cos(angles(m)))/...
                          (impedence_layers(m+1)*cos(angles(m))+impedence_layers(m)*cos(angles(m+1)));
end

%--------------------------------------------------------------------------
% Matrix equations for the parallel and perpendicular components
% Below is the initial matrix for the first interface in the structure
eqn_P = [1/transmit_coeff_P(1) reflect_coeff_P(1)/transmit_coeff_P(1);...
          reflect_coeff_P(1)/transmit_coeff_P(1) 1/transmit_coeff_P(1)];
eqn_S  = [1/transmit_coeff_S(1) reflect_coeff_S(1)/transmit_coeff_S(1);...
          reflect_coeff_S(1)/transmit_coeff_S(1) 1/transmit_coeff_S(1)];
      
%If there are more interfaces than just 1, the following code is used to
%combine the rest of the interfaces and propigation components, as two
%interfaces are required for any "bounces" and changes in the propigation.
max_interfaces = length(transmit_coeff_P);
if(max_interfaces > 1)
    for m = 2:(max_interfaces)
       propagation_S = [exp(1i*k_z_layers(m)*layer_thickness(m)) 0 ;...
                        0 exp(-1i*k_z_layers(m)*layer_thickness(m))];
       propagation_P= [exp(1i*k_z_layers(m)*layer_thickness(m)) 0 ;...
                        0 exp(-1i*k_z_layers(m)*layer_thickness(m))];
       interface_S   = [1/transmit_coeff_S(m) reflect_coeff_S(m)/transmit_coeff_S(m);...
                        reflect_coeff_S(m)/transmit_coeff_S(m) 1/transmit_coeff_S(m)];
       interface_P  = [1/transmit_coeff_P(m) reflect_coeff_P(m)/transmit_coeff_P(m);...
                        reflect_coeff_P(m)/transmit_coeff_P(m) 1/transmit_coeff_P(m)];
       eqn_S = eqn_S*propagation_S*interface_S;
       eqn_P= eqn_P*propagation_P*interface_P;
    end
end
eqn_S = eqn_S*[Eto_T;0];
eqn_P = eqn_P*[Eto_II;0];

%Solve the final matricies

Eto_perpendicular = double(Eio/(eqn_S(1)/Eto_T));
Ero_perpendicular = double((eqn_S(2)/Eto_T)*Eto_perpendicular);

Eto_parallel = double(Eio/(eqn_P(1)/Eto_II));
Ero_parallel = double((eqn_P(2)/Eto_II)*Eto_parallel);

%Add the solved reflection and transmission coefficients to the overall
%data matrix to get a graph of the coefficients as the incident wavelength
%is varied.
ref_data_reflection_S(wavelength) = abs(Ero_perpendicular)^2;
ref_data_transmission_S(wavelength) = abs(Eto_perpendicular)^2;
ref_data_reflection_P(wavelength) = abs(Ero_parallel)^2;
ref_data_transmission_P(wavelength) = abs(Eto_parallel)^2;
totals_S(wavelength) = abs(Ero_perpendicular)^2 + abs(Eto_perpendicular)^2;
totals_P(wavelength) = abs(Ero_parallel)^2 + abs(Eto_parallel)^2;
end

%Plot the data on two seperate graphs.
figure()
subplot(2, 1, 1)
plot(range, ref_data_reflection_S)
title("Perpendicular Reflection Magnitude")
%xlabel("Incident wavelength (meters)")
ylabel("Perpendicular Reflection Magnitude")
subplot(2, 1, 2)
plot(range, ref_data_transmission_S)
title("Perpendicular Transmission Magnitude")
%xlabel("Incident wavelength (meters)")
ylabel("Perpendicular Transmission Magnitude")

figure()
subplot(2, 1, 1)
plot(range, ref_data_reflection_P)
title("Parallel Reflection Magnitude")
%xlabel("Incident wavelength (meters)")
ylabel("Parallel Reflection Magnitude")
subplot(2, 1, 2)
plot(range, ref_data_transmission_P)
title("Parallel Transmission Magnitude")
%xlabel("Incident wavelength (meters)")
ylabel("Parallel Transmission Magnitude")

figure()
subplot(2, 1, 1)
plot(range, totals_S)
title("Sum of Transmission and Reflection for the Perpendicular Case")
subplot(2, 1, 2)
plot(range, totals_P)
title("Sum of Transmission and Reflection for the Parallel Case")