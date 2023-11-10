% Incorporation Model for continuous Villermaux-Dushman reaction
% determination of micro mixing time

% Author: Torben Frey  torben.frey@tuhh.de
% sources:
% original model [Fou96], doi:10.1016/S0009-2509(96)00340-5                      
% modified model [Ari21], doi:10.1016/j.cherd.2021.09.010                           

%% assumptions 
% reaction system
% (I)    TRIS  +  H+           -->  TRISH+          buffer reaction
% (II)   5I-  +  IO3-  +  6H+  -->  3I2  +  3H2O    dushman reaction
% (III)  I2  +  I-            <-->  I3-             symproportionation
% species indices
% (1) = H+      (2) = TRIS          (3) = TRISH+
% (4) = I-      (5) = IO3-          (6) = I2 
% (7) = H2O     (8) = I3-

%% preamble

% USER INPUT: incorporation function
fcn = 'exp';        % 'lin' or 'exp' incorporation function
mdl = 'arian';      % 'fournier' (original) or 'arian' (modified)

% USER INPUT: volume flows V1 and V2
V(1) = 2;     	    % mL/min buffer volume flow 
V(2) = V(1) ;	    % mL/min acid volume flow

% USER INPUT: start values of iterative detemination of micromixing time
t_start = 0.2;      % t in s: initial guess value for micro mixing time

% USER INPUT: experimental measurements
Xs_exp = 0.0546;    % segregation index (experiment) 
I3_exp = 0.1374E-03;    % triiodide concentration (experiment) [mol/L]

% USER INPUT: initial concentration of each species
c0(1) = 0.03;       % mol/L (H+)
c0(2) = 0.0898;     % mol/L (TRIS) 
c0(3) = 0.0898;     % mol/L (TRISH+)
c0(4) = 0.03197;    % mol/L (I-) 
c0(5) = 6.34E-3;    % mol/L (IO3-) 
c0(6) = 0.0;        % mol/L (I2)
c0(7) = 0.0;        % mol/L (H2O) / n.a.
c0(8) = 0.0;        % mol/L (I3-)

% USER INPUT: Setting tolerances for ODE solver
optionsODE = odeset('RelTol',1E-20,'AbsTol',1E-21); % ODE tolerances
optionsFMIN = optimset('TolFun',1E-20); % fminsearch tolerances

%% determination of function minima

warning('off','all')

tic
% reaction time tr1 is found when Xs_exp and Xs (model) coincide
tr1 = fminsearch(@(tm)(minimum_xs(tm,c0,V,Xs_exp,fcn,mdl,optionsODE)),t_start,optionsFMIN);
% micro mixing time tm1 is found when I3_exp and I3 (model) coincide
tm1 = fminsearch(@(tm)(minimum_I3(tm,c0,V,I3_exp,fcn,mdl,optionsODE)),t_start,optionsFMIN); 
% obtain time arrays of the model for tm1
[Xs1,n1,v1,t1] = VD_incorporation_model(tm1,c0,V,Xs_exp,fcn,mdl,optionsODE); 
toc