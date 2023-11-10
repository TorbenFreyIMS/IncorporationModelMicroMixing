% Incorporation model for continuous Villermaux-Dushman reaction
% Determining Delta of model and experiment

% Author: Torben Frey  torben.frey@tuhh.de
% sources:
% original model [Fou96], doi:10.1016/S0009-2509(96)00340-5                       
% modified model [Ari21], doi:10.1016/j.cherd.2021.09.010              

% FUNCTION INPUTS: 
% tm: initial guess for micro mixing time [s]
% c0: initial concentration array nxm = 1x8 [mol/L]
% V: volume flow array V(1):V1, V(2):V2 [mL/h]
% Xs_exp: experimentally measured segregation index
% fcn: incorporation function setting. options: "lin" or "exp"
% mdl: incorporation model setting. options: "fournier" or "arian"
% options: ODE15s solver options

% OUTPUT:
% deltaXs: squared difference of experimental and model segregation index

function [deltaXs] = minimum_xs(tm,c0,V,Xs_exp,fcn,mdl,options)
% convert flow rate to L/s 
V = V/1000/60;     	    % L/s buffer volume flow 
% initial molar flow for each species n=c*V=[mol/L*L/s]
n0(1) = c0(1) * V(2);   % mol/s (H+)
n0(2) = c0(2) * V(1);   % mol/s (TRIS)
n0(3) = c0(3) * V(1);   % mol/s (TRISH+)
n0(4) = c0(4) * V(1);   % mol/s (I-)
n0(5) = c0(5) * V(1);   % mol/s (IO3-)
n0(6) = c0(6) * V(1);   % mol/s (I2)
n0(7) = c0(7) * V(1);   % mol/s (H20) / n.a.
n0(8) = c0(8) * V(1);   % mol/s (I3-)
% set initial conditions
nt = [n0(1) 0 0 0 0 0 0 0]; % mol/s
% set integration interval [0 tend] depending on model
if strcmp(mdl,"fournier") 
    if strcmp(fcn,"lin")
        % tend for linear function (Fournier)
        tend = V(1) / V(2) * tm;
    elseif strcmp(fcn,"exp") % fcn == "exp"
        % tend for exponential function (Fournier)
        tend = log( (V(1) + V(2)) / V(2) ) * tm;    
    else
        disp('Invalid format of "fcn". Options: "lin" or "exp". Proceed with linear incorporation function.')
        % fallback: tend for linear function (Fournier)
        tend = V(1) / V(2) * tm;
    end
elseif strcmp(mdl,"arian") 
    if strcmp(fcn,"lin")
        % tend for linear function (Arian)
        tend = tm;                                  
    elseif strcmp(fcn,"exp") 
        % tend for exponential function (Arian)
        tend = 5 * tm;                              
    else
        disp('Invalid format of "fcn". Options: "lin" or "exp". Proceed with linear incorporation function.')
        % fallback: tend for linear function (Arian)
        tend = tm; 
    end
else
    disp('Invalid format of "mdl". Options: "fournier" or "arian". Proceed with modified incorporation model.')
    if strcmp(fcn,"lin") % fcn == "lin"
        % fallback: tend for linear function (Arian)
        tend = tm;
    elseif strcmp(fcn,"exp") % fcn == "exp"
        % fallback: tend for exponential function (Arian)
        tend = 5 * tm;
    else
        disp('Invalid format of "fcn". Options: "lin" or "exp". Proceed with linear incorporation function.')
        % fallback: tend for linear function (Arian)
        tend = tm; 
    end
end

% ODE solver
[~,n] = ode15s(@ODE_solver_frey, [0 tend], nt, options, n0, V, tm, fcn, mdl);


% stoichiometric yield Yst of the Dushman reaction
% Yst = 6*IO3(0) / ( 6*IO3(0) + TRIS(0) )
Yst = 6*n0(5) / (6*n0(5) + n0(2));
% actual yield Y of the dushman reaction
% Y = 2*(I2(end)+I3(end)) / H(0)
Y =  2*(n(end,6)+n(end,8)) / n0(1);
% segregation index from incorporation model
Xs =  Y / Yst;

deltaXs = (Xs-Xs_exp)^2;

fprintf('(Xs-Xs_exp)%c = %.3g at tr = %f ms\n',178,deltaXs,tm*1000)

end