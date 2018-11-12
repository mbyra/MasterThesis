% This file is part of script calculating QFI(t) and state coefficients of 
% given state in cavity with displaced mirror.
%
% Main script. Sets all parameters, loads initial states' declarations,
% creates necessary directories, and runs computation for given ranges of
% time and g. Displays and exports QFI(t) after simulation for given g;
% displays and exports multiplot QFI(t) for all gs.
%
% Author: Marcin Byra, UW
% email: marcin.byra1@gmail.com
% 09/2018

global omegaM omega0 f maxSteps accuracy initialNbar g t chartsVisibility...
    N debug
% load states declarations
load states.mat


% *************************************************************************


% CONTROL OF THE SIMULATION

% Parameters - the same or close to those in Mathematica in order to 
% compare results
omegaM = 0.3;
omega0 = 10;
f = 0.001;
maxSteps = 100;
accuracy = 10;

% Select state to generate QFI(t) and coefficients(t) for.
initialState = state5Constr;

% Set on/off to view/hide figs with coefficient charts (every 0.25 t unit)
chartsVisibility = 'on';

% Toggle visibility of debug info
debug = true;


% *************************************************************************


% ACTUAL CODE

N = length(initialState);
initialNbar = (N-1)/2;
format shortG;

% create directory for storing figures and jpegs:
fn = fullfile('figures');
if ~exist('f', 'dir')
   warning('Creating directory figures/'); 
   mkdir(fn);
end
    
time = 1:0.05:3; 
g_list = [0.01, 0.04, 0.08, 0.15, 0.3, 0.5]; % coupling constant
qfi_values = zeros(length(g_list), length(time));


fprintf(sprintf('Starting a simulation of %-dimensional state:\n', N));
display(initialState);
fprintf('accuracy = %d, max steps = %d, omega0 = %.2f, omegaM = %.2f\n',...
    accuracy, maxSteps, omega0, omegaM);
% main loop generating evoultion for various values of g
for a = 6:length(g_list)
    g = g_list(a);
    fprintf('\n g = %.2f...\n', g);
    
    for b = 1:length(time) % set plot density; units 2Pi/omega_0
        t = time(b);
        fprintf('\t t = %.2f...\n', t);
        
        [newState, qfi_values(a, b)] = evolution(initialState, 'initial state');
        
%         TODO
%         evolution(fliplr(result), 'reversed initial state');
%         evolution(initialState/2 + fliplr(result)/2,...
%             'superposition of initial and reversed result state')
    end
    
    % show QFI(t) plot and export png to subdirectory 'figures'
    fh = figure('Name', sprintf('qfi(t) for N = %d and g = %d', N, g*100));
    title(sprintf('QFI(t)for %d-dimensional initial state, constrained', N));
    semilogy(time, qfi_values(a, :));
    xlabel('t [2\pi/{\omega}m]');
    ylabel('QFI');
    legend(sprintf("g = %.2f", g), 'Location', 'southeast');
    saveFmt = '\\figures\\qfi_vs_t_dim_%d_g_%03d.png';
    saveas(fh,[pwd sprintf(saveFmt, N, g*100)]);
end

% display QFI(t) multiplot for every g and expert png to subdir 'figures'
fh = figure('Name', sprintf('qfi(t) for N = %d and all g', N));
semilogy(time, qfi_values(1, :), time, qfi_values(2, :), time, ...
    qfi_values(3, :), time, qfi_values(4, :), time, qfi_values(5, :),...
    time, qfi_values(6, :));
title(sprintf('QFI(t)for %d-dimensional initial state, constrained', N));
xlabel('t [2\pi/{\omega}m]');
ylabel('QFI');
legend('g = 0.01', 'g = 0.04','g = 0.08', 'g = 0.15', 'g = 0.3', ...
    'g = 0.5', 'Location', 'southeast');
saveDirFmt = '\\figures\\multiplot_qfi_vs_t_dim_%d_g_%03d.png';
saveas(fh,[pwd sprintf(saveDirFmt, N, g*100)]);




% Simulates evolution of initial state for a given t and parameters. Prints
% info to command window and checks if calculation was valid; if not,
% changed qfi to NaN to avoid plotting this value.
% Returns a pair: final state and its qfi.
function [state, qfi] = evolution(currentState, stateStr) 
    global omegaM omega0 f maxSteps accuracy initialNbar g t
    
    % Writing because it can be initial, reversed or superposition state
    fprintf(sprintf('\t\tRunning evolution of %s...', stateStr));
    
    [state, qfi, steps] = calculateOptimalQFI(currentState);
    fprintf(' Done.\n');

    if steps > 0
        fprintf('\t\tAccuracy reached after %d steps.\n', steps);
    else
        fprintf('\t\tCould not reach accuracy after %d steps.\n', maxSteps);
    end
    fprintf('\t\tFound QFI = %f.\n\n', real(qfi));

    %just to avoid creating ugly plots when calculation was invalid:
    if steps <= 0 || real(qfi) < 10
        qfi = NaN;
    end

    qfi = real(qfi);
end
