% This file is part of script calculating QFI(t) and state coefficients of 
% given state in cavity with displaced mirror.
%
% Declarations of initial states used to generate time dependencies and
% coefficient charts. The same states were used in appropiate plots in
% different implementation (Mathematica).
%
% This script has to be run every time it is modified to export updated
% variables to states.mat file.
%
% Author: Marcin Byra, UW
% email: marcin.byra1@gmail.com
% 10/2018

% The same initial states as in Mathematica (previous, unconstrained)
state5 = [0.2, 0.1i, 0.3i, 0.1, sqrt(1 - 0.2^2 - 0.3^2 - 2*0.1^2)];
state5 = state5/norm(state5);

state10 = [0.2, 0.1i, 0.3i, 0.1, 0.5, 0.2, 0.1i, 0.7i, 0.5, 0.05i];
state10 = state10/norm(state10);

state15 = [0.2, 0.1i, 0.3i, 0.1, 0.5, 0.2, 0.1i, 0.7i, 0.5, 0.05i, 0.8, ...
    0.1, 0.05i, 0.4i, 0.2];
state15 = state15/norm(state15);

% The same initial states as in Mathematica (current, fulfilling constraints)
state5Constr = [-0.17143-0.483256i, -0.184004-0.466691i,...
    0.0000238363-0.0000768779i, 0.341161+0.277375i, -0.525566-0.125976i];
state15Constr = [0.265834-0.0252019i, -0.00523809-0.00495137i,...
                 0.294972-0.260994i, -0.0907457+0.293633i,...
                 0.0140262-0.0181147i, -0.0833171+0.212157i,...
                 0.0121167-0.00771775i, 0.258398+0.324393i,...
                 0.3176+0.0554487i, 0.00118436-0.00185454i,...
                 0.420994-0.0446258i, -0.00185298+0.00507704i,...
                 0.0198546-0.193464i, -0.0628111+0.0426242i,...
                 0.124764 +0.3347i];

% Run this script after every modification.
save('states.mat', 'state5', 'state10', 'state15', 'state5Constr', ...
     'state15Constr');