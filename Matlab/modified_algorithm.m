omegaM = 0.3;
omega0 = 10;
f = 0.001;
maxSteps = 100;
accuracy = 10;

state5 = [0.2, 0.1i, 0.3i, 0.1, sqrt(1 - 0.2^2 - 0.3^2 - 2*0.1^2)];
state5 = state5/norm(state5);
state10 = [0.2, 0.1i, 0.3i, 0.1, 0.5, 0.2, 0.1i, 0.7i, 0.5, 0.05i];
state101 = state10/norm(state10);
state15 = [0.2, 0.1i, 0.3i, 0.1, 0.5, 0.2, 0.1i, 0.7i, 0.5, 0.05i, 0.8, ...
    0.1, 0.05i, 0.4i, 0.2];
state15 = state15/norm(state15);

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

current_state = state5Constr;

N = length(current_state);
initialNbar = (N-1)/2;
fprintf("N = %d\n", N);
format shortG;

% create directory for storing figures and jpegs:
fn = fullfile('figures');
if ~exist('f', 'dir')
   warning('Creating directory figures/'); 
   mkdir(fn);
end

time = 1:0.05:3;
g_list = [0.01, 0.04, 0.08, 0.15, 0.3, 0.5];
qfi_values = zeros(length(g_list), length(time));
% 
% fig = figure;
% title('QFI(t) for 5-dimensional initial state, constrained');
for a = 5:length(g_list)
    g = g_list(a);
    fprintf('Starting loop for g = %.2f...\n', g);
    for b = 1:length(time)
        t = time(b);
        fprintf('\tRunning calculations for t = %.2f...', t);
        [rho, qfi, steps] = calculateOptimalQFI(f, g, (t + 0.0000001)*2*pi/omegaM,...
            omegaM, omega0, maxSteps, accuracy, current_state, initialNbar, t);
        fprintf('Done.\n');
        if steps > 0
            fprintf('\t\tAccuracy reached after %d steps.\n', steps);
        else
            fprintf('\t\tCould not reach accuracy after %d steps.\n', maxSteps);
%             qfi = NaN;
        end
        fprintf('\t\tFound QFI = %f.\n.', real(qfi));
        qfi_values(a, b) = real(qfi);
        
%         fprintf('\tRunning calculations for t = %.2f... REVERSED STATE', t);
%         [rho, qfi, steps] = calculateOptimalQFI(f, g, (t + 0.0000001)*2*pi/omegaM,...
%             omegaM, omega0, maxSteps, accuracy, fliplr(current_state), initialNbar, t);
%         fprintf('Done.\n');
%         if steps > 0
%             fprintf('\t\tAccuracy reached after %d steps.\n', steps);
%         else
%             fprintf('\t\tCould not reach accuracy after %d steps.\n', maxSteps);
% %             qfi = NaN;
%         end
%         fprintf('\t\tFound QFI = %f.\n.', real(qfi));
%         
%         fprintf('\tRunning calculations for t = %.2f... SUPERPOSITION OF RESULTS', t);
%         [rho, qfi, steps] = calculateOptimalQFI(f, g, (t + 0.0000001)*2*pi/omegaM,...
%             omegaM, omega0, maxSteps, accuracy, fliplr(current_state), initialNbar, t);
%         fprintf('Done.\n');
%         if steps > 0
%             fprintf('\t\tAccuracy reached after %d steps.\n', steps);
%         else
%             fprintf('\t\tCould not reach accuracy after %d steps.\n', maxSteps);
% %             qfi = NaN;
%         end
%         fprintf('\t\tFound QFI = %f.\n.', real(qfi));
        
    end
    fig = figure(a);
%     subplot(2,3,a);
    semilogy(time, qfi_values(a, :));
    title(sprintf('dim=%d, constrained, iterations:%d', N, maxSteps));
    xlabel('t [2\pi/{\omega}m]'); % x-axis label
    ylabel('QFI');
    legend(sprintf("g = %.2f", g), 'Location', 'southeast');
    saveas(fig,[pwd sprintf('\\figures\\qfi_vs_t_dim_%d_g_%03d_iterations_%d.png'...
    , N, g*100, maxSteps)]);
end

saveas(fig,[pwd sprintf('\\figures\\multisubplot_qfi_vs_t_dim_%d_g_%03d_iterations_%d.png'...
    , N, g*100, maxSteps)]);

fig = figure;
semilogy(time, qfi_values(1, :), time, qfi_values(2, :), time, ...
    qfi_values(3, :), time, qfi_values(4, :), time, qfi_values(5, :),...
    time, qfi_values(6, :));
title(sprintf('dim=%d, constrained, iterations:%d', N, maxSteps));
xlabel('t [2\pi/{\omega}m]'); % x-axis label
ylabel('QFI');
legend('g = 0.01', 'g = 0.04','g = 0.08', 'g = 0.15', 'g = 0.3', ...
    'g = 0.5', 'Location', 'southeast');
saveas(fig,[pwd sprintf('\\figures\\multiplot_qfi_vs_t_dim_%d_g_%03d_iterations_%d.png'...
    , N, g*100, maxSteps)]);


