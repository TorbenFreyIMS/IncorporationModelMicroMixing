% Incorporation model for continuous Villermaux-Dushman reaction
% Generation of temporal arrays of Xs, n, v

% Author: Torben Frey (-4124), torben.frey@tuhh.de
% sources:  original incorporation model (Fournier, DOI: 10.1016/S0009-2509(96)00340-5)
%           modified incorporation model (Arian, DOI: 10.1016/j.cherd.2021.09.010) 

% INPUTS: 
% tm: known micro mixing time [s]
% c0: initial concentration array nxm = 1x8 [mol/L]
% V: volume flow array V(1):V1, V(2):V2 [mL/h]
% NOT USED (~) 
% fcn: incorporation function setting. options: "lin" or "exp"
% mdl: incorporation model setting. options: "fournier" or "arian"
% options: ODE15s solver options

% OUTPUTS:
% Xs: Segregation index array over time/length [-]
% n: species flux array array over time/length [mol/s]
% v: incorporation volume flow V(2) array over time/length [L/s]
% t: time/length array [s]

function [Xs, n, v, t] = VD_incorporation_model(tm, c0, V, ~, fcn, mdl, options)

% SI units of volume flow
V = V/1000/60;     	        % L/s buffer volume flow 
% initial molar flow for each species n=c*V=[mol/L*L/s]
n0(1) = c0(1) * V(2);   % mol/s (H+)
n0(2) = c0(2) * V(1);   % mol/s (TRIS)
n0(3) = c0(3) * V(1);   % mol/s (TRISH+)
n0(4) = c0(4) * V(1);   % mmol/s (I-)
n0(5) = c0(5) * V(1);   % mol/s (IO3-)
n0(6) = c0(6) * V(1);   % mol/s (I2)
n0(7) = c0(7) * V(1);   % mol/s (H20) / n.a.
n0(8) = c0(8) * V(1);   % mol/s (I3-)

% Initial conditions
nt = [n0(1) 0 0 0 0 0 0 0]; % mol/s

% Integration interval [0 tend]
if strcmp(mdl,"fournier") % mdl == "fournier"
    if strcmp(fcn,"lin") % fcn == "lin"
        tend = V(1) / V(2) * tm;                    % tend for linear function (Fournier et al.)
    elseif strcmp(fcn,"exp") % fcn == "exp"
        tend = log( (V(1) + V(2)) / V(2) ) * tm;    % tend for exponential function (Fournier et al.)
    else
        warning('Invalid format of "fcn". Options: "lin" or "exp". Proceed with linear incorporation function.')
        tend = V(1) / V(2) * tm;
    end
elseif strcmp(mdl,"arian") % mdl == "arian"
    if strcmp(fcn,"lin") % fcn == "lin"
        tend = tm;                                  % tend for linear function (Arian 2021)
    elseif strcmp(fcn,"exp") % fcn == "exp"
        tend = 5 * tm;                              % tend for exponential function (Arian 2021)
    else
        warning('Invalid format of "fcn". Options: "lin" or "exp". Proceed with linear incorporation function.')
        tend = tm; 
    end
else
    warning('Invalid format of "mdl". Options: "fournier" or "arian". Proceed with modified incorporation model.')
    if strcmp(fcn,"lin") % fcn == "lin"
        tend = tm;                                  % tend for linear function (Arian 2021)
    elseif strcmp(fcn,"exp") % fcn == "exp"
        tend = 5 * tm;                              % tend for exponential function (Arian 2021)
    else
        warning('Invalid format of "fcn". Options: "lin" or "exp". Proceed with linear incorporation function.')
        tend = tm; 
    end
end

% ODE solver
[t,n] = ode15s(@ODE_solver_frey, [0 tend], nt, options, n0, V, tm, fcn, mdl);

% maximum yield Yst of the Dushman reaction (equivalent to maximum segregation)
Yst = 6*n0(5) / (6*n0(5) + n0(2));
% actual yield Y of the dushman reaction
Y =  2*(n(:,6)+n(:,8)) / n0(1);
% segregation index from incorporation model
Xs =  Y / Yst;

if strcmp(mdl,"fournier") % mdl == "fournier"
    if strcmp(fcn,"lin") % fcn == "lin"
        v = V(2) * (1+t/tm(end)); 
    elseif strcmp(fcn,"exp") % fcn == "exp"
        v = V(2) * exp(t/tm(end));
    else
        warning('Invalid format of "fcn". Options: "lin" or "exp". Proceed with linear incorporation function.')
        v = V(2) * (1+t/tm(end)); 
    end
elseif strcmp(mdl,"arian") % mdl == "arian"
    if strcmp(fcn,"lin") % fcn == "lin"
        v = V(2) + V(1) * (t/tm(end));
    elseif strcmp(fcn,"exp") % fcn == "exp"
        v = V(2) + V(1) * (1-exp(-t/tm(end)));
    else
        warning('Invalid format of "fcn". Options: "lin" or "exp". Proceed with linear incorporation function.')
        v = V(2) + V(1) * (t/tm(end));
    end
else
    warning('Invalid format of "mdl". Options: "fournier" or "arian". Proceed with modified incorporation model.')
    if strcmp(fcn,"lin") % fcn == "lin"
        v = V(2) + V(1) * (t/tm(end));
    elseif strcmp(fcn,"exp") % fcn == "exp"
        v = V(2) + V(1) * (1-exp(-t/tm(end)));
    else
        warning('Invalid format of "fcn". Options: "lin" or "exp". Proceed with linear incorporation function.')
        v = V(2) + V(1) * (t/tm(end));
    end
end

end