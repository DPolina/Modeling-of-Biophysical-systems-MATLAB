clc;
clear all;
close all;

% Начальные условия
T0 = 0; % Начальное время
T = 20; % Конечное время
N = 10000; % Количество симуляций
m = 50; % Количество интервалов для разбиения временного отрезка
dt = T / m; % Шаг разбиения

% f2 = (2.5)
f1 = 1;
prob = 0.005;
lambda_m = -log(prob) * f1;
disp(['lambda_m = ', num2str(lambda_m)]);

param_1 = f1^2;
param_2 = 1 / 2;
sigma = (f1^2 / 4)^(1/4);
alpha = (1 / 2) * param_2;

% Генерация потока Кокса (стационарный)

flow = cell(N, 1); % Поток

for i = 1:N

    t = []; % Реализация стац. потока Пуассона

    t_next = T0 - log(rand) / lambda_m;
    while(t_next >= T0 && t_next <= T) 
        t = [t t_next];
        t_prev = t_next;
        t_next = t_prev-log(rand)/lambda_m;
    end

    % Независимые идентичные стационарные нормальные случайные процессы 
    x1 = zeros(1,length(t)); % реализация одного из случайных процессов
    x2 = zeros(1,length(t));

    fi1 = randn(1, length(t)); % независимые реализации нормальной случайной величины 
    fi2 = randn(1, length(t));
    x1(1) = sigma * fi1(1); % (2.6) фи - норм. распред. [0, 1]
    x2(1) = sigma * fi2(1);

    % (2.6)
    for k = 2:length(t)
        x1(k) = x1(k-1)*exp(-alpha*(t(k)-t(k-1))) + sigma*fi1(k)*sqrt(1-exp(-2*alpha*(t(k)-t(k-1))));
        x2(k) = x2(k-1)*exp(-alpha*(t(k)-t(k-1))) + sigma*fi2(k)*sqrt(1-exp(-2*alpha*(t(k)-t(k-1))));
    end
    
    % (2.4) 
    xi = x1.^2 + x2.^2;

    % Прореживание
    % Для каждого значения экспоненциального случайного процесса xi моделируем 
    % реализацию случайной величины, равномерно распределенной на интервале
    % [0; lambda_m].
    y = lambda_m * rand([1, length(t)]); 
    t_fin = t(y < xi);
    
    %disp('t_fin')
    %disp(t_fin);

    flow{i} = t_fin;
end

% Оценка интенсивности 

X = T0 + dt/2 : dt : T-dt/2; % центры карманов

hist_vect = zeros(1, m);

for i = 1:N
    hist_temp = hist(flow{i}, X);
    hist_vect = hist_vect + hist_temp;
end

% Вычисляем среднее количество событий
empirical = hist_vect / (N * dt);
teoretical = f1;

figure();
hold on;
plot(X, empirical, 'r', 'DisplayName', 'Экспериментальная интенсивность')
plot(X, teoretical * ones(1, length(X)), 'b', 'DisplayName', 'Теоретическая интенсивность')
title("Оценка интенсивности ");
xlabel("T");
ylabel("Интенсивность");
legend()
hold off;

% Оценка корреляционной функции

% Инициализация матрицы для корреляционных функций (Кол-во реализаций потока N (строки) х
% кол-во карманов (столбцы)). Эта матрица будет использоваться для хранения 
% рассчитанных значений корреляционной функции для каждой реализации и каждого временного сдвига.
K = zeros(1, m);

% (1.8) - для стац. процессов.
for L = 1:length(flow)
    k = hist(flow{L}, X);
    for i = 1:m
        summ = 0;
        if i == 1
            for j = 1:m
                summ = summ + (k(j) * (k(j) - 1))/(dt^2);
            end
            K(i) = K(i) + summ/m;
        else
            for j = 1:m - (i-1)
                summ = summ + (k(j)*k(j+i-1))/(dt^2);
            end
            K(i) = K(i) + summ/(m - (i-1));
        end
    end
end

K = K / length(flow);
K_theor = param_1 * (1 + exp(-X * param_2)); % на основании (2.1)
% Моментные функции потока Кокса тождественно совпадают с моментными 
% функциями исходного случайного процесса 𝜉(𝑡). 

figure();
hold on;
plot(X, K, 'r', 'DisplayName', 'Экспериментальная корреляция');
plot(X, K_theor, 'b', 'DisplayName', 'Теоретическая корреляция');
title('График корреляционной функции');
legend();
hold off;
title("Оценка корреляционной функции");

% Распределение количества событий
% Вер-ть того, что в реализации потока будет n событий.

tc = 2;
A_p = tc / 10; % Интервал времени для подсчета событий.

A_temp = cell(N, 1); % хранение событий в каждом интервале времени для каждой симуляции

for n = 1:length(A_temp) 
    % Выбираются события, произошедшие до A_p.
    A_temp{n} = flow{n}(flow{n} < A_p);
end

% Нахождение максимального числа событий
max_events = length(A_temp{1});
for i = 2:length(A_temp)
    % Определяем максимальное количество событий в любом интервале времени 
    % среди всех симуляций.
    if length(A_temp{i}) > max_events
        max_events = length(A_temp{i});
    end
end

% Вычисление эмпирического распределения событий
P = zeros(1, max_events + 1);
% Для каждого возможного числа событий от 0 до max_num:
for i = 0:max_events
    % Подсчитывается количество симуляций с таким числом событий,
    % и вычисляется относительная частота этого числа событий.
    for n = 1:length(A_temp)
        if length(A_temp{n}) == i
            P(i+1) = P(i+1) + 1;
        end
    end
    P(i+1) = P(i+1) / length(A_temp);
end

% Вычисление теоретического распределения событий
P_theor = zeros(1, max_events+1);
m = f1 * A_p;
for k = 0:max_events
    P_theor(k+1) = m^k / ((m + 1) ^ (k+1)); % (2.3)
end

figure();
hold on;
plot(0:1:max_events, P, 'r', 'DisplayName', 'Экспериментальное распределение');
plot(0:1:max_events, P_theor, 'b', 'DisplayName', 'Теоретическое распределение');
legend();
hold off;
title("Распределение количества событий");
