clc 
clear
close all

% Load the error data
error8 = load("total_error_8x8.mat");
error8 = error8.total_error;
error16 = load("total_error_16x16.mat");
error16 = error16.total_error;
error32 = load("total_error_32x32.mat");
error32 = error32.total_error;
error_list = [error8 error16 error32];

% Mesh sizes
h = [1/8 1/16 1/32];

% Calculate slopes
slopes = diff(log(error_list)) ./ diff(log(h));

% Log-log plot

figure;
fig = gcf;
loglog(h, error_list, '-o', 'LineWidth', 1.5);
hold on;

% Annotate points and slopes
for i = 1:length(h)
    
    % Show the slope on the intervals
    if i < length(h)
        % Midpoint for slope annotation
        mid_h = sqrt(h(i) * h(i+1));
        mid_error = sqrt(error_list(i) * error_list(i+1));
        text(mid_h, mid_error, sprintf('Slope = %.2f', slopes(i)), 'HorizontalAlignment', 'center');
    end
end

% Enhance plot aesthetics
grid on;
xlabel('Mesh Size (h)', 'FontSize', 10);
ylabel('Log Error', 'FontSize', 10);
title('Log-Log Plot of Error vs. Mesh Size', 'FontSize', 10);
set(gca, 'XScale', 'log', 'YScale', 'log');
exportgraphics(fig,"error_plot.png", "Resolution",300)

