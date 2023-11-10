% ODE15s solver for incorporation model for continuous Villermaux-Dushman

% Author: Torben Frey, torben.frey@tuhh.de
% sources:  
% original model [Fou96], doi:10.1016/S0009-2509(96)00340-5           
% modified model [Ari21], doi:10.1016/j.cherd.2021.09.010                
% acid-base reaction constants [Owe34], doi:10.1021/ja01323a014              
% triiodide equilibrium constants [Rua86], doi:10.1021/j100409a034          

%% FUNCTION INPUTS:

% t: start/end time interval [t1 t2]
% nt: concentration array start conditions nxm = 1x8
% n0: concentration array initial conditions nxm = 1x8
% V: volume array V(1):V1, V(2):V2
% tm: micro mixing time 
% fcn: incorporation function setting. options: "lin" or "exp"
% mdl: incorporation model setting. options: "fournier" or "arian"

function [dndt] =ODE_solver_frey(t, nt, n0, V, tm, fcn, mdl)

%% definition of incorporation laws
if mdl == "fournier"
    if fcn == "lin"    
        % linear incorporation function (Fournier)
        g = 1+t/tm;
        dgdt = 1/tm;
    elseif fcn == "exp"    
        % exponential incorporation function (Fournier et al.)
        g = exp(t/tm);     
        dgdt = exp(t/tm)/tm;
    else
        disp('Invalid format of "fcn". Options: "lin" or "exp". Proceed with linear incorporation function.')
        % fallback: linear incorporation function (Fournier)
        g = 1+t/tm;             
        dgdt = 1/tm;
    end
    % volume growth of V2(t) (Fournier)
    v = V(2) * g;
elseif mdl == "arian"
    if fcn == "lin"   
        % linear incorporation function (Arian)
        g = t/tm;           
        dgdt = 1/tm;
    elseif fcn == "exp"      
        % exponential incorporation function (Arian)
        g = 1-exp(-t/tm);  
        dgdt = exp(-t/tm)/tm;   
    else
        disp('Invalid format of "fcn". Options: "lin" or "exp". Proceed with linear incorporation function.')
        % fallback: linear incorporation function (Arian)
        g = t/tm;
        dgdt = 1/tm;
    end
    % volume growth of V2(t) (Arian)
    v = V(2) + V(1) * g;    
else
    disp('Invalid format of "mdl". Options: "fournier" or "arian". Proceed with modified incorporation model.')
    if fcn == "lin"         
        % fallback: linear incorporation function (Arian)
        g = t/tm;           
        dgdt = 1/tm;
    elseif fcn == "exp"         
        % fallback: exponential incorporation function (Arian)
        g = 1-exp(-t/tm);
        dgdt = exp(-t/tm)/tm;   
    else
        disp('Invalid format of "fcn". Options: "lin" or "exp". Proceed with linear incorporation function.')
        % fallback: linear incorporation function (Arian)
        g = t/tm;
        dgdt = 1/tm;
    end
    % fallback: volume growth of V2(t) (Arian)
    v = V(2) + V(1) * g;    
end


%% rate constants & reaction rates

% reaction system
% (I)    TRIS  +  H+          <-->  TRISH+          buffer reaction
% (II)   5I-  +  IO3-  +  6H+  -->  3I2  +  3H2O    Dushman reaction
% (III)  I2  +  I-            <-->  I3-             triiodide equilibrium
% species indices
% (1) = H+      (2) = TRIS          (3) = TRISH+
% (4) = I-      (5) = IO3-          (6) = I2 
% (7) = H2O     (8) = I3-

% time-dependant ionic strength for perchloric acid and TRIS buffer
% Z = 1/2 *sum(charge^2*n) / V2
% overall charge = time dependent molar flux nt PLUS inlet molar flux n0 of counter-ions (e.g. ClO4-, K+)
Z = 0.5/v * (n0(1)+nt(1)+nt(3)+n0(4)+nt(4)+n0(5)+nt(5)+nt(8));
      
% reaction constants k(i)
% k(1) acid-buffer reaction (Owen 1934) in L/mol/s
k(1) = 1e11;                                        
% k(2) acid-buffer backwards reaction (Owen 1934) in s^-1
k(2) = 8.71e2;                                       
% k(3) Dushman reaction (Arian 2021) in L^4/mol^4/s
k(3) = 1.37E9*10^(-1.93*sqrt(Z)/(1+sqrt(Z))+0.40*Z);
% k(4) triiode equilibrium L/mol/s (Ruasse 1986) in L/mol/s
k(4) = 5.6e9;   
% k(5) triiode equilibrium backwards L/mol/s (Ruasse 1986) in s^-1
k(5) = 7.5e6;                                        

% acid-buffer reaction rate (I)
r1 = k(1)/v * nt(1) * nt(2) - k(2) * nt(3);   
 % Dushman reaction rate (II)
r2 = k(3)/v^4 * nt(1)^2 * nt(4)^2 * nt(5);   
% triiodide equilibrium reaction reate (III)
r3 = k(4)/v * nt(4) * nt(6) - k(5) * nt(8);   

%% time dependent species arrays (ODE)

if mdl == "fournier"         
% balance equations of original incorporation model (Fournier)
    dndt = zeros(length(nt),1); % preallocate dndt
    dndt(1) = -r1 - 6*r2;                             % mol/s (H+)
    dndt(2) = -r1             + n0(2)*dgdt*V(2)/V(1); % mol/s (TRIS)
    dndt(3) =  r1             + n0(3)*dgdt*V(2)/V(1); % mol/s (TRISH+)
    dndt(4) =     - 5*r2 - r3 + n0(4)*dgdt*V(2)/V(1); % mol/s (I-)
    dndt(5) =     - r2        + n0(5)*dgdt*V(2)/V(1); % mol/s (IO3-)
    dndt(6) =       3*r2 - r3;                        % mol/s (I2)
    dndt(7) =       3*r2;                             % mol/s (H20)
    dndt(8) =              r3;                        % mol/s (I3-)
else                                        
% balance equations of modified incorporation model (Arian)
    dndt = zeros(length(nt),1); % preallocate dndt
    dndt(1) = -r1 - 6*r2;                   % mol/s (H+)
    dndt(2) = -r1             + n0(2)*dgdt; % mol/s (TRIS)
    dndt(3) =  r1             + n0(3)*dgdt; % mol/s (TRISH+)
    dndt(4) =     - 5*r2 - r3 + n0(4)*dgdt; % mol/s (I-)
    dndt(5) =     - r2        + n0(5)*dgdt; % mol/s (IO3-)
    dndt(6) =       3*r2 - r3;              % mol/s (I2)
    dndt(7) =       3*r2;                   % mol/s (H20)
    dndt(8) =              r3;              % mol/s (I3-)
end

return