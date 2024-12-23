clc
close all
clear all

% Load distances
load('distances1.mat', 'rn');
n = length(rn);

% Define block size
block_size = 100;
num_rows = n - block_size;

% Create shift table
PT = zeros(num_rows, block_size);

for i = 1:num_rows
    PT(i, :) = rn(i:i+block_size-1);
end

train_size = floor(0.8 * num_rows);
test_size = num_rows - train_size;

% Split data
P = PT(1:train_size, :);
T = PT(train_size+1:end, :);

% Create network with two hidden layers of 5 neurons.
hidden = [5 5];
net = newff(P,T,hidden);

% Simulate network and its output plotted against the targets before 
% training.

Y = sim(net,P);
figure()
hold on;
plot(T(:, 1), '-', 'DisplayName', 'T');
plot(Y(:, 1), 'o', 'DisplayName', 'Y');
legend();

%xlabel('Input');
%ylabel('Output');
%legend('Target', 'Network Output');
title('Before training');

% Network is trained for 50 epochs. 
net.trainParam.epochs = 20;
net = train(net,P,T);
Y = sim(net,P);

figure()
hold on;
plot(T(:, 1), '-', 'DisplayName', 'T');
plot(Y(:, 1), 'o', 'DisplayName', 'Y');
legend();
title('After 50 epochs');

% Testing
X = PT(test_size+1:end, :);
Yg = sim(net,X);

figure()
hold on;
plot(X(:, 1), '-', 'DisplayName', 'X');
plot((length(X(:, 1)) + 1):(length(X(:, 1)) + length(Yg(:, 1))), Yg(:, 1), '--', 'DisplayName', 'Yg');
legend();
title('Testing');

whos;