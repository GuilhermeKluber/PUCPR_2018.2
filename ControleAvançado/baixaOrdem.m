clear all; %Limpa as variaveis salvas no Workspace
close all; %Fecha os gr�ficos abertos na tela
clc;       %Limpa o console do Matlab
s = tf('s');  %Configura a letra 's' como uma variavels complexa


%parametros
tempo_inicial=0;
tempo_amostragem = 0.01;
tempo_final = 20;
G = 1 / (s+1)^8;  %Fun��o de Tranferencia


%Tempos: Inicial, Amostragem e Final
t = tempo_inicial:tempo_amostragem:tempo_final;  % Define o tempo inicial (0), o passo (0.01) e o tempo final (20)

%Degrau
U = heaviside(t); %Cria um vetor degrau, s� � usado para plotar no grafico

%Grafico: Resposta
step(G,t);  %Mostra a resposta ao degrau aplicado na fun��o de transferencia G, durante o tempo t
hold on;   %Comando para incluir mais inform��oes no gr�fico

%Grafico: Degrau
plot(t,U, 'r');  %plota o degrau criado em U, somente para visualiza��o
hold on;



%Valores da resposta ao degrau
[y] = step(G,t); %Cria um vetor com os valores da resposta ao degrau aplicado na tf G

y_min = min(y); %valor inicial da Resposta
y_max = max(y);%valor final da Resposta

u_min = 0; %menor valor do degrau
u_max = max(U); %Valor final do degrau


%reta tangente
y1 = 0.48*y_max;
y2 = 0.52*y_max;
x1 = t(1,find(y>y1,1));
x2 = t(1,find(y>y2,1));

syms x 


m = (y2 - y1) / (x2 - x1);
f(x) = m*(x - x1) +y1;


%Grafico: Reta Tangente
ezplot(f, [tempo_inicial,tempo_final]);


title('Fun��o degrau'); %define o titulo do gr�fico
xlabel('tempo')
ylabel('resposta')
legend('Respota','Degrau','Reta Tangente')
axis([(tempo_inicial-1) tempo_final (y_min-1) (y_max+1)]);  %configura os eixos que devem ser exibidos no gr�fico [x(min) x(max) y(min) y(max)]


%Inc�gnitas
K =0;
tau =0;
teta = 0;

%Modelo da fun��o de primeiro grau a ser encontrada pelos diferentes
%m�todos

H = exp(-teta*s)*(K/(tau*s + 1));


%Definindo o valor de K, todos os m�todos usam a mesma f�rmula para o ganho
%K

%M�todo de Ziegler-Nichols (1942)
%M�todo de H�gglund (1991)
%M�todo de Smith (1985)
%M�todo de Sundaresan e Krishnaswamy (1977)
%M�todo de Nishikawa (1984)

%K = varia��o da saida dividida pela varia��o da entrada

K = (y_max- y_min)/(u_max-u_min);



%Definindo o valor de tau e teta
%M�todo de Ziegler-Nichols (1942)
%M�todo de H�gglund (1991)

eqn = f(x)==0;
cross_bottom =  double (solve(eqn,x));
plot(cross_bottom,0,'ro');


eqn = f(x)==y_max;
cross_top=  double (solve(eqn,x));
plot(cross_top,y_max,'ro');

zig_tau = cross_top - cross_bottom;
zig_teta = cross_bottom-tempo_inicial;

hag_tau =  t(1,find(y>(0.632*y_max),1))  - cross_bottom;
hag_teta = zig_teta;


%M�todo de Smith (1985)
smith_tau = 1.5*(t(1,find(y>(0.632*y_max),1)) - t(1,find(y>(0.283*y_max),1)));
smith_teta = t(1,find(y>(0.632*y_max),1))-smith_tau;

%M�todo de Sundaresan e Krishnaswamy (1977)
sund_tau = 0.67*(t(1,find(y>(0.853*y_max),1)) - t(1,find(y>(0.353*y_max),1)));
sund_teta = 1.3*(t(1,find(y>(0.353*y_max),1)))-0.29*(t(1,find(y>(0.853*y_max),1)));

%M�todo de Nishikawa (1984)
%Integral da diferen�a entre o valor da saido no infinito e o valor da
%saido no tempo.
A0 = trapz(t,(y_max-y(:,1))); %faz a integral num�rica da diferen�a entre as curvas, obtendo a area A0
t0 = A0/y_max;
new_t = tempo_inicial:tempo_amostragem:t0;
t0_size = numel(new_t);
A1 = trapz(new_t,(y(1:t0_size,1)));

nish_tau = A1/0.368*y_max;
nish_teta = t0 - nish_tau;

%Gr�fico Ziegler
figure('Name','Ziegler');
step(G,t);  %Mostra a resposta ao degrau aplicado na fun��o de transferencia G, durante o tempo t
hold on;   %Comando para incluir mais inform��oes no gr�fico
tau = zig_tau;
teta = zig_teta;
H = exp(-teta*s)*(K/(tau*s + 1));
step(H,t);

real = step(G,t);
obtida = step(H,t);

% Mean Squared Error (MSE) � objetivo � um menor MSE
MSE = sum((real-obtida).^2) / (length(real));
desempenho = ['Ziegler  -  MSE=', num2str(MSE),'  K=',num2str(K),'  tau=',num2str(tau),'  teta=',num2str(teta)];

title(desempenho); %define o titulo do gr�fico
xlabel('tempo')
ylabel('resposta')
legend('Original','Ziegler')
axis([(tempo_inicial-1) tempo_final (y_min-1) (y_max+1)]);  %configura os eixos que devem ser exibidos no gr�fico [x(min) x(max) y(min) y(max)]


%Gr�fico H�gglund
figure('Name','H�gglund');
step(G,t);  %Mostra a resposta ao degrau aplicado na fun��o de transferencia G, durante o tempo t
hold on;   %Comando para incluir mais inform��oes no gr�fico
tau = hag_tau;
teta = hag_teta;
H = exp(-teta*s)*(K/(tau*s + 1));
step(H,t);

obtida = step(H,t);

% Mean Squared Error (MSE) � objetivo � um menor MSE
MSE = sum((real-obtida).^2) / (length(real));
desempenho = ['H�gglund  -  MSE=', num2str(MSE),'  K=',num2str(K),'  tau=',num2str(tau),'  teta=',num2str(teta)];

title(desempenho); %define o titulo do gr�fico

xlabel('tempo')
ylabel('resposta')
legend('Original','H�gglund')
axis([(tempo_inicial-1) tempo_final (y_min-1) (y_max+1)]);  %configura os eixos que devem ser exibidos no gr�fico [x(min) x(max) y(min) y(max)]


%Gr�fico Smith
figure('Name','Smith');
step(G,t);  %Mostra a resposta ao degrau aplicado na fun��o de transferencia G, durante o tempo t
hold on;   %Comando para incluir mais inform��oes no gr�fico
tau = smith_tau;
teta = smith_teta;
H = exp(-teta*s)*(K/(tau*s + 1));
step(H,t);

obtida = step(H,t);

% Mean Squared Error (MSE) � objetivo � um menor MSE
MSE = sum((real-obtida).^2) / (length(real));
desempenho = ['Smith  -  MSE=', num2str(MSE),'  K=',num2str(K),'  tau=',num2str(tau),'  teta=',num2str(teta)];

title(desempenho); %define o titulo do gr�fico

xlabel('tempo')
ylabel('resposta')
legend('Original','Smith')
axis([(tempo_inicial-1) tempo_final (y_min-1) (y_max+1)]);  %configura os eixos que devem ser exibidos no gr�fico [x(min) x(max) y(min) y(max)]


%Gr�fico Sundaresan
figure('Name','Sundaresan');
step(G,t);  %Mostra a resposta ao degrau aplicado na fun��o de transferencia G, durante o tempo t
hold on;   %Comando para incluir mais inform��oes no gr�fico
tau = sund_tau;
teta = sund_teta;
H = exp(-teta*s)*(K/(tau*s + 1));
step(H,t);

obtida = step(H,t);

% Mean Squared Error (MSE) � objetivo � um menor MSE
MSE = sum((real-obtida).^2) / (length(real));
desempenho = ['Sundaresan  -  MSE=', num2str(MSE),'  K=',num2str(K),'  tau=',num2str(tau),'  teta=',num2str(teta)];

title(desempenho); %define o titulo do gr�fico

xlabel('tempo')
ylabel('resposta')
legend('Original','Sundaresan')
axis([(tempo_inicial-1) tempo_final (y_min-1) (y_max+1)]);  %configura os eixos que devem ser exibidos no gr�fico [x(min) x(max) y(min) y(max)]



%Gr�fico Nishikawa
figure('Name','Nishikawa');
step(G,t);  %Mostra a resposta ao degrau aplicado na fun��o de transferencia G, durante o tempo t
hold on;   %Comando para incluir mais inform��oes no gr�fico
tau = nish_tau;
teta = nish_teta;
H = exp(-teta*s)*(K/(tau*s + 1));
step(H,t);

obtida = step(H,t);

% Mean Squared Error (MSE) � objetivo � um menor MSE
MSE = sum((real-obtida).^2) / (length(real));
desempenho = ['Nishikawa  -  MSE=', num2str(MSE),'  K=',num2str(K),'  tau=',num2str(tau),'  teta=',num2str(teta)];

title(desempenho); %define o titulo do gr�fico

xlabel('tempo')
ylabel('resposta')
legend('Original','Nishikawa')
axis([(tempo_inicial-1) tempo_final (y_min-1) (y_max+1)]);  %configura os eixos que devem ser exibidos no gr�fico [x(min) x(max) y(min) y(max)]






