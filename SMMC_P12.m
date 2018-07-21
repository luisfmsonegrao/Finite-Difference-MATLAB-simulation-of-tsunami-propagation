%Nota: organizei o trabalho um pouco mal, pois em vez de escrever funções
%reutilizáveis construí um bloco de código enorme. a leitura dos
%comentários ajuda a compreender a lógica do código
%(por exemplo, é nos comentários que refiro o desaparecimento dos termos
%viscosos - as derivadas cruzadas de v - e de compressibilidade - a
%divergência de v)
R=10;
g = 9.81;
A = 3*pi()/4;
N = 100;
dx = 2*R/N;
dy = 2*R/N;
x = linspace(0, 2*R, N);
y = linspace(0, 2*R, N);
[X,Y] = meshgrid(x,y);
dt = 0.002;
B = zeros(N,N);
dB_dx = B;
dB_dy = B;
%% fundo
for i = 1:N
    for j = 1:N
        r = sqrt( (j*dx-R)^2 + (i*dy-R)^2 );
        if ( r <= R/3 )
            B(i,j) = - 10;
            dB_dx(i,j) = 0;
            dB_dy(i,j) = 0;
        elseif r <= R
            B(i,j) = -10 * (cos((A) * (r-R/3)/R))^2;
            dB_dx(i,j) = 20*((j*dx-R)/r)*(A/R)*sin(A*(r-R/3)/R)*cos(A*(r-R/3)/R);%derivadas analiticas de B em ordem a x e y
            dB_dy(i,j) = 20*((i*dy-R)/r)*(A/R)*sin(A*(r-R/3)/R)*cos(A*(r-R/3)/R);%derivadas analiticas de B em ordem a x e y
        else 
            B(i,j)=0;
            dB_dx(i,j) = 0;
            dB_dy(i,j) = 0;
        end
    end
end
dB_dx(:,1) = 0; 
dB_dy(1,:) = 0;
figure(01);
g1=surf(X, Y, B);
set(g1,'LineStyle','none');
title('Seabottom topology');
xlabel('x(m)');
ylabel('y(m)');
%% Altura inicial
h = exp (- ((X-R).^2 + (Y-R).^2));%condicao inicial
figure(02);
h1 = surf(X,Y,h);
set(h1,'LineStyle','none');
title('Initial Wave Height');
xlabel('x(m)');
ylabel('y(m)');
vx = zeros(N,N);
vy = zeros(N,N);
vxnew = vx;%vou guardar as atualizacoes de vx, vy e h numa nova matriz e so depois de iterar sobre todos os pontos atualizar a matriz antiga
vynew = vy;
hnew = h;
ti = 0;
tmax = 1.35; %% a partir daqui a onda reflete e começa lentamente a divergir;
hmax=zeros(tmax/dt);% par guardar a altura maxima da onda em funcao de t
t_hmax = zeros(tmax/dt);
vol = zeros(tmax/dt); %para guardar a estimativa do volume deslocado em funcao de t
t_count = 1;
graph_count = 0;
while (ti < tmax)
    
    for i = 1:N
        
        for j = 1:N % usei condicoes fronteira periodicas, p. exemplo se i = 1 em vez de v(i-1,j) uso v(N,j)
            
            if ((j*dx-R)^2 + (i*dx-R)^2 < R^2) %fora desta regiao o fluido passa a ser inviscido e incompressivel, pelo que desaparecem as derivadas cruzadas da velocidade no calculo de vx e vy e a divergencia da velocidade no calculo de h!!!
                if i == 1 %exclui j's afastados do centro (j= 1,2,N-1 e N)
                    
                    vxnew(i,j) = vx(i,j) - dt*(vx(i,j)*(-vx(i,j+2)+8*vx(i,j+1)-8*vx(i,j-1)+vx(i,j-2))/(12*dx) + vy(i,j)*(-vx(i+2,j)+8*vx(i+1,j)-8*vx(N,j)+vx(N-1,j))/(12*dy) + g*((-h(i,j+2)+8*h(i,j+1)-8*h(i,j-1)+h(i,j-2))/(12*dx)));
                    vynew(i,j) = vy(i,j) - dt*(vy(i,j)*(-vy(i+2,j)+8*vy(i+1,j)-8*vy(N,j)+vy(N-1,j))/(12*dy) + vx(i,j)*(-vy(i,j+2)+8*vy(i,j+1)-8*vy(i,j-1)+vy(i,j-2))/(12*dx) + g*((-h(i+2,j)+8*h(i+1,j)-8*h(N,j)+h(N-1,j))/(12*dx)));
                    if h(i,j) > B(i,j) %condicao para garantir que o nivel de agua nao e mais baixa que o fundo oceanico
                        hnew(i,j) = h(i,j) - dt*(((-vx(i,j+2)+8*vx(i,j+1)-8*vx(i,j-1)+vx(i,j-2)) + (-vy(i+2,j)+8*vy(i+1,j)-8*vy(N,j)+vy(N-1,j)))*(h(i,j)-B(i,j))/(12*dy) + vx(i,j)*((h(i,j+1)-h(i,j-1))/(2*dx)-dB_dx(i,j)) + vy(i,j)*((h(i+1,j)-h(N,j))/(2*dx)-dB_dy(i,j)));
                    end
                elseif i == 2
                    vxnew(i,j) = vx(i,j) - dt*(vx(i,j)*(-vx(i,j+2)+8*vx(i,j+1)-8*vx(i,j-1)+vx(i,j-2))/(12*dx) + vy(i,j)*(-vx(i+2,j)+8*vx(i+1,j)-8*vx(i-1,j)+vx(N,j))/(12*dy) + g*((-h(i,j+2)+8*h(i,j+1)-8*h(i,j-1)+h(i,j-2))/(12*dx)));
                    vynew(i,j) = vy(i,j) - dt*(vy(i,j)*(-vy(i+2,j)+8*vy(i+1,j)-8*vy(i-1,j)+vy(N,j))/(12*dy) + vx(i,j)*(-vy(i,j+2)+8*vy(i,j+1)-8*vy(i,j-1)+vy(i,j-2))/(12*dx) + g*((-h(i+2,j)+8*h(i+1,j)-8*h(i-1,j)+h(N,j))/(12*dx)));
                    if h(i,j) > B(i,j)
                        hnew(i,j) = h(i,j) - dt*(((-vx(i,j+2)+8*vx(i,j+1)-8*vx(i,j-1)+vx(i,j-2)) + (-vy(i+2,j)+8*vy(i+1,j)-8*vy(i-1,j)+vy(N,j)))*(h(i,j)-B(i,j))/(12*dy) + vx(i,j)*((h(i,j+1)-h(i,j-1))/(2*dx)-dB_dx(i,j)) + vy(i,j)*((h(i+1,j)-h(i-1,j))/(2*dx)-dB_dy(i,j)));
                    end
                elseif i == N-1
                    vxnew(i,j) = vx(i,j) - dt*(vx(i,j)*(-vx(i,j+2)+8*vx(i,j+1)-8*vx(i,j-1)+vx(i,j-2))/(12*dx) + vy(i,j)*(-vx(1,j)+8*vx(i+1,j)-8*vx(i-1,j)+vx(i-2,j))/(12*dy) + g*((-h(i,j+2)+8*h(i,j+1)-8*h(i,j-1)+h(i,j-2))/(12*dx)));
                    vynew(i,j) = vy(i,j) - dt*(vy(i,j)*(-vy(1,j)+8*vy(i+1,j)-8*vy(i-1,j)+vy(i-2,j))/(12*dy) + vx(i,j)*(-vy(i,j+2)+8*vy(i,j+1)-8*vy(i,j-1)+vy(i,j-2))/(12*dx) + g*((-h(1,j)+8*h(i+1,j)-8*h(i-1,j)+h(i-2,j))/(12*dx)));
                    if h(i,j) > B(i,j)
                        hnew(i,j) = h(i,j) - dt*(((-vx(i,j+2)+8*vx(i,j+1)-8*vx(i,j-1)+vx(i,j-2)) + (-vy(1,j)+8*vy(i+1,j)-8*vy(i-1,j)+vy(i-2,j)))*(h(i,j)-B(i,j))/(12*dy) + vx(i,j)*((h(i,j+1)-h(i,j-1))/(2*dx)-dB_dx(i,j)) + vy(i,j)*((h(i+1,j)-h(i-1,j))/(2*dx)-dB_dy(i,j)));
                    end
                elseif i == N
                    vxnew(i,j) = vx(i,j) - dt*(vx(i,j)*(-vx(i,j+2)+8*vx(i,j+1)-8*vx(i,j-1)+vx(i,j-2))/(12*dx) + vy(i,j)*(-vx(2,j)+8*vx(1,j)-8*vx(i-1,j)+vx(i+2,j))/(12*dy) + g*((-h(i,j+2)+8*h(i,j+1)-8*h(i,j-1)+h(i,j-2))/(12*dx)));
                    vynew(i,j) = vy(i,j) - dt*(vy(i,j)*(-vy(2,j)+8*vy(1,j)-8*vy(i-1,j)+vy(i-2,j))/(12*dy) + vx(i,j)*(-vy(i,j+2)+8*vy(i,j+1)-8*vy(i,j-1)+vy(i,j-2))/(12*dx) + g*((-h(2,j)+8*h(1,j)-8*h(i-1,j)+h(i-2,j))/(12*dx)));
                    if h(i,j) > B(i,j)
                         hnew(i,j) = h(i,j) - dt*(((-vx(i,j+2)+8*vx(i,j+1)-8*vx(i,j-1)+vx(i,j-2)) + (-vy(2,j)+8*vy(1,j)-8*vy(j-1,j)+vy(j-2,j)))*(h(i,j)-B(i,j))/(12*dy) + vx(i,j)*((h(i,j+1)-h(i,j-1))/(2*dx)-dB_dx(i,j)) + vy(i,j)*((h(1,j)-h(i-1,j))/(2*dx)-dB_dy(i,j)));
                    end
                elseif j == 1
                    vxnew(i,j) = vx(i,j) - dt*(vx(i,j)*(-vx(i,j+2)+8*vx(i,j+1)-8*vx(i,N)+vx(i,N-1))/(12*dx) + vy(i,j)*(-vx(i+2,j)+8*vx(i+1,j)-8*vx(i-1,j)+vx(i-2,j))/(12*dy) + g*((-h(i,j+2)+8*h(i,j+1)-8*h(i,N)+h(i,N-1))/(12*dx)));
                    vynew(i,j) = vy(i,j) - dt*(vy(i,j)*(-vy(i+2,j)+8*vy(i+1,j)-8*vy(i-1,j)+vy(i-2,j))/(12*dy) + vx(i,j)*(-vy(i,j+2)+8*vy(i,j+1)-8*vy(i,N)+vy(i,N-1))/(12*dx) + g*((-h(i+2,j)+8*h(i+1,j)-8*h(i-1,j)+h(i-2,j))/(12*dx)));
                    if h(i,j) > B(i,j)
                         hnew(i,j) = h(i,j) - dt*(((-vx(i,j+2)+8*vx(i,j+1)-8*vx(i,N)+vx(i,N-1)) + (-vy(i+2,j)+8*vy(i+1,j)-8*vy(i-1,j)+vy(i-2,j)))*(h(i,j)-B(i,j))/(12*dy) + vx(i,j)*((h(i,j+1)-h(i,N))/(2*dx)-dB_dx(i,j)) + vy(i,j)*((h(i+1,j)-h(i-1,j))/(2*dx)-dB_dy(i,j)));
                    end
                elseif j == 2
                    vxnew(i,j) = vx(i,j) - dt*(vx(i,j)*(-vx(i,j+2)+8*vx(i,j+1)-8*vx(i,j-1)+vx(i,N))/(12*dx) + vy(i,j)*(-vx(i+2,j)+8*vx(i+1,j)-8*vx(i-1,j)+vx(i-2,j))/(12*dy) + g*((-h(i,j+2)+8*h(i,j+1)-8*h(i,j-1)+h(i,N))/(12*dx)));
                    vynew(i,j) = vy(i,j) - dt*(vy(i,j)*(-vy(i+2,j)+8*vy(i+1,j)-8*vy(i-1,j)+vy(i-2,j))/(12*dy) + vx(i,j)*(-vy(i,j+2)+8*vy(i,j+1)-8*vy(i,j-1)+vy(i,N))/(12*dx) + g*((-h(i+2,j)+8*h(i+1,j)-8*h(i-1,j)+h(i-2,j))/(12*dy)));
                    if h(i,j) > B(i,j)
                        hnew(i,j) = h(i,j) - dt*((((-vx(i,j+2)+8*vx(i,j+1)-8*vx(i,j-1)+vx(i,N)) + (-vy(i+2,j)+8*vy(i+1,j)-8*vy(i-1,j)+vy(i-2,j)))*(h(i,j)-B(i,j))/(12*dy) + vx(i,j)*((h(i,j+1)-h(i,j-1))/(2*dx) - dB_dx(i,j)) + vy(i,j)*((h(i+1,j)-h(i-1,j))/(2*dy)-dB_dy(i,j))));
                    end
                elseif j == N-1
                    vxnew(i,j) = vx(i,j) - dt*(vx(i,j)*(-vx(i,1)+8*vx(i,j+1)-8*vx(i,j-1)+vx(i,j-2))/(12*dx) + vy(i,j)*(-vx(i+2,j)+8*vx(i+1,j)-8*vx(i-1,j)+vx(i-2,j))/(12*dy) + g*((-h(i,1)+8*h(i,j+1)-8*h(i,j-1)+h(i,j-2))/(12*dx)));
                    vynew(i,j) = vy(i,j) - dt*(vy(i,j)*(-vy(i+2,j)+8*vy(i+1,j)-8*vy(i-1,j)+vy(i-2,j))/(12*dy) + vx(i,j)*(-vy(i,1)+8*vy(i,j+1)-8*vy(i,j-1)+vy(i,j-2))/(12*dx) + g*((-h(i+2,j)+8*h(i+1,j)-8*h(i-1,j)+h(i-2,j))/(12*dx)));
                    if h(i,j) > B(i,j)
                        hnew(i,j) = h(i,j) - dt*(((-vx(i,1)+8*vx(i,j+1)-8*vx(i,j-1)+vx(i,j-2)) + (-vy(i+2,j)+8*vy(i+1,j)-8*vy(i-1,j)+vy(i-2,j)))*(h(i,j)-B(i,j))/(12*dy) + vx(i,j)*((h(i,j+1)-h(i,j-1))/(2*dx) - dB_dx(i,j)) + vy(i,j)*((h(i+1,j)-h(i-1,j))/(2*dy)-dB_dy(i,j)));
                    end
                elseif j == N
                    vxnew(i,j) = vx(i,j) - dt*(vx(i,j)*(-vx(i,2)+8*vx(i,1)-8*vx(i,j-1)+vx(i,j-2))/(12*dx)+vy(i,j)*(-vx(i+2,j)+8*vx(i+1,j)-8*vx(i-1,j)+vx(i-2,j))/(12*dy) + g*((-h(i,2)+8*h(i,1)-8*h(i,j-1)+h(i,j-2)))/(12*dx));
                    vynew(i,j) = vy(i,j) - dt*(vy(i,j)*(-vy(i+2,j)+8*vy(i+1,j)-8*vy(i-1,j)+vy(i-2,j))/(12*dy) + vx(i,j)*(-vy(i,2)+8*vy(i,1)-8*vy(i,j-1)+vy(i,j-2))/(12*dx) + g*((-h(i+2,j)+8*h(i+1,j)-8*h(i-1,j)+h(i-2,j))/(12*dy)));
                    if h(i,j) > B(i,j)
                         hnew(i,j) = h(i,j) - dt*(((-vx(i,2)+8*vx(i,1)-8*vx(i,j-1)+vx(i,j-2)) + (-vy(i+2,j)+8*vy(i+1,j)-8*vy(i-1,j)+vy(i-2,j)))*(h(i,j)-B(i,j))/(12*dy) + vx(i,j)*((h(i,1)-h(i,j-1))/(2*dx)-dB_dx(i,j)) + vy(i,j)*((h(i+1,j)-h(i-1,j))/(2*dx)-dB_dy(i,j)));
                    end
                else
                    vxnew(i,j) = vx(i,j) - dt*(vx(i,j)*(-vx(i,j+2)+8*vx(i,j+1)-8*vx(i,j-1)+vx(i,j-2))/(12*dx)+vy(i,j)*(-vx(i+2,j)+8*vx(i+1,j)-8*vx(i-1,j)+vx(i-2,j))/(12*dy) + g*((-h(i,j+2)+8*h(i,j+1)-8*h(i,j-1)+h(i,j-2))/(12*dx)));
                    vynew(i,j) = vy(i,j) - dt*(vy(i,j)*(-vy(i+2,j)+8*vy(i+1,j)-8*vy(i-1,j)+vy(i-2,j))/(12*dy)+vx(i,j)*(-vy(i,j+2)+8*vy(i,j+1)-8*vy(i,j-1)+vy(i,j-2))/(12*dx) + g*((-h(i+2,j)+8*h(i+1,j)-8*h(i-1,j)+h(i-2,j))/(12*dy)));
                    if h(i,j) > B(i,j)
                       hnew(i,j) = h(i,j) - dt*(((-vx(i,j+2)+8*vx(i,j+1)-8*vx(i,j-1)+vx(i,j-2)) + (-vy(i+2,j)+8*vy(i+1,j)-8*vy(i-1,j)+vy(i-2,j)))*(h(i,j)-B(i,j))/(12*dy) + vx(i,j)*((h(i,j+1)-h(i,j-1))/(2*dx) - dB_dx(i,j)) + vy(i,j)*((h(i+1,j)-h(i-1,j))/(2*dy)-dB_dy(i,j)));
                    end
                end
            else % fora de R < = 10
                if i == 1
                    if j == 1
                        vxnew(i,j) = vx(i,j) - dt*(vx(i,j)*(-vx(i,j+2)+8*vx(i,j+1)-8*vx(i,N)+vx(i,N-1))/(12*dx) + g*((-h(i,j+2)+8*h(i,j+1)-8*h(i,N)+h(i,N-1))/(12*dx)));
                        vynew(i,j) = vy(i,j) - dt*(vy(i,j)*(-vy(i+2,j)+8*vy(i+1,j)-8*vy(N,j)+vy(N-1,j))/(12*dy) + g*((-h(i+2,j)+8*h(i+1,j)-8*h(N,j)+h(N-1,j))/(12*dy)));
                        if h(i,j) > B(i,j)
                             hnew(i,j) = h(i,j) - dt*(vx(i,j)*((h(i,j+1)-h(i,N))/(2*dx)-dB_dx(i,j)) + vy(i,j)*((h(i+1,j)-h(N,j))/(2*dx)-dB_dy(i,j)));
                        end
                    elseif j == 2
                        vxnew(i,j) = vx(i,j) - dt*(vx(i,j)*(-vx(i,j+2)+8*vx(i,j+1)-8*vx(i,j-1)+vx(i,N))/(12*dx) + g*((-h(i,j+2)+8*h(i,j+1)-8*h(i,j-1)+h(i,N))/(12*dx)));
                        vynew(i,j) = vy(i,j) - dt*(vy(i,j)*(-vy(i+2,j)+8*vy(i+1,j)-8*vy(N,j)+vy(N-1,j))/(12*dy) + g*((-h(i+2,j)+8*h(i+1,j)-8*h(N,j)+h(N-1,j))/(12*dy)));
                        if h(i,j) > B(i,j)
                            hnew(i,j) = h(i,j) - dt*(vx(i,j)*((h(i,j+1)-h(i,j-1))/(2*dx)-dB_dx(i,j)) + vy(i,j)*((h(i+1,j)-h(N,j))/(2*dx)-dB_dy(i,j)));
                        end
                    elseif j == N-1
                        vxnew(i,j) = vx(i,j) - dt*(vx(i,j)*(-vx(i,1)+8*vx(i,j+1)-8*vx(i,j-1)+vx(i,j-2))/(12*dx) + g*((-h(i,1)+8*h(i,j+1)-8*h(i,j-1)+h(i,j-2))/(12*dx)));
                        vynew(i,j) = vy(i,j) - dt*(vy(i,j)*(-vy(i+2,j)+8*vy(i+1,j)-8*vy(N,j)+vy(N-1,j))/(12*dy) + g*((-h(i+2,j)+8*h(i+1,j)-8*h(N,j)+h(N-1,j))/(12*dy)));
                        if h(i,j) > B(i,j)
                             hnew(i,j) = h(i,j) - dt*(vx(i,j)*((h(i,j+1)-h(i,j-1))/(2*dx)-dB_dx(i,j)) + vy(i,j)*((h(i+1,j)-h(N,j))/(2*dx)-dB_dy(i,j)));
                        end
                    elseif j == N
                        vxnew(i,j) = vx(i,j) - dt*(vx(i,j)*(-vx(i,2)+8*vx(i,1)-8*vx(i,j-1)+vx(i,j-2))/(12*dx) + g*((-h(i,2)+8*h(i,1)-8*h(i,j-1)+h(i,j-2))/(12*dx)));
                        vynew(i,j) = vy(i,j) - dt*(vy(i,j)*(-vy(i+2,j)+8*vy(i+1,j)-8*vy(N,j)+vy(N-1,j))/(12*dy) + g*((-h(i+2,j)+8*h(i+1,j)-8*h(N,j)+h(N-1,j))/(12*dy)));
                        if h(i,j) > B(i,j)
                             hnew(i,j) = h(i,j) - dt*(vx(i,j)*((h(i,1)-h(i,j-1))/(2*dx)-dB_dx(i,j)) + vy(i,j)*((h(i+1,j)-h(N,j))/(2*dx)-dB_dy(i,j)));
                        end
                    end
                elseif i == 2
                    if j == 1
                        vxnew(i,j) = vx(i,j) - dt*(vx(i,j)*(-vx(i,j+2)+8*vx(i,j+1)-8*vx(i,N)+vx(i,N-1))/(12*dx) + g*((-h(i,j+2)+8*h(i,j+1)-8*h(i,N)+h(i,N-1))/(12*dx)));
                        vynew(i,j) = vy(i,j) - dt*(vy(i,j)*(-vy(i+2,j)+8*vy(i+1,j)-8*vy(i-1,j)+vy(N-1,j))/(12*dy) + g*((-h(i+2,j)+8*h(i+1,j)-8*h(i-1,j)+h(N,j))/(12*dy)));
                        if h(i,j) > B(i,j)
                             hnew(i,j) = h(i,j) - dt*(vx(i,j)*((h(i,j+1)-h(i,N))/(2*dx)-dB_dx(i,j)) + vy(i,j)*((h(i+1,j)-h(i-1,j))/(2*dx)-dB_dy(i,j)));
                        end
                    elseif j == 2
                        vxnew(i,j) = vx(i,j) - dt*(vx(i,j)*(-vx(i,j+2)+8*vx(i,j+1)-8*vx(i,j-1)+vx(i,N))/(12*dx) + g*((-h(i,j+2)+8*h(i,j+1)-8*h(i,j-1)+h(i,N))/(12*dx)));
                        vynew(i,j) = vy(i,j) - dt*(vy(i,j)*(-vy(i+2,j)+8*vy(i+1,j)-8*vy(i-1,j)+vy(N,j))/(12*dy) + g*((-h(i+2,j)+8*h(i+1,j)-8*h(i-1,j)+h(N,j))/(12*dy)));
                        if h(i,j) > B(i,j)
                            hnew(i,j) = h(i,j) - dt*(vx(i,j)*((h(i,j+1)-h(i,j-1))/(2*dx)-dB_dx(i,j)) + vy(i,j)*((h(i+1,j)-h(i-1,j))/(2*dx)-dB_dy(i,j)));
                        end
                    elseif j == N-1
                        vxnew(i,j) = vx(i,j) - dt*(vx(i,j)*(-vx(i,1)+8*vx(i,j+1)-8*vx(i,j-1)+vx(i,j-2))/(12*dx) + g*((-h(i,1)+8*h(i,j+1)-8*h(i,j-1)+h(i,j-2))/(12*dx)));
                        vynew(i,j) = vy(i,j) - dt*(vy(i,j)*(-vy(i+2,j)+8*vy(i+1,j)-8*vy(i-1,j)+vy(N,j))/(12*dy) + g*((-h(i+2,j)+8*h(i+1,j)-8*h(i-1,j)+h(N,j))/(12*dy)));
                        if h(i,j) > B(i,j)
                            hnew(i,j) = h(i,j) - dt*(vx(i,j)*((h(i,j+1)-h(i,j-1))/(2*dx)-dB_dx(i,j)) + vy(i,j)*((h(i+1,j)-h(i-1,j))/(2*dx)-dB_dy(i,j)));
                        end
                    elseif j == N
                        vxnew(i,j) = vx(i,j) - dt*(vx(i,j)*(-vx(i,2)+8*vx(i,1)-8*vx(i,j-1)+vx(i,j-2))/(12*dx) + g*((-h(i,2)+8*h(i,1)-8*h(i,j-1)+h(i,j-2))/(12*dx)));
                        vynew(i,j) = vy(i,j) - dt*(vy(i,j)*(-vy(i+2,j)+8*vy(i+1,j)-8*vy(i-1,j)+vy(N,j))/(12*dy) + g*((-h(i+2,j)+8*h(i+1,j)-8*h(i-1,j)+h(N,j))/(12*dy)));
                        if h(i,j) > B(i,j)
                            hnew(i,j) = h(i,j) - dt*(vx(i,j)*((h(i,1)-h(i,j-1))/(2*dx)-dB_dx(i,j)) + vy(i,j)*((h(i+1,j)-h(i-1,j))/(2*dx)-dB_dy(i,j)));
                        end
                    end
                elseif i == N-1
                    if j == 1
                        vxnew(i,j) = vx(i,j) - dt*(vx(i,j)*(-vx(i,j+2)+8*vx(i,j+1)-8*vx(i,N)+vx(i,N-1))/(12*dx) + g*((-h(i,j+2)+8*h(i,j+1)-8*h(i,N)+h(i,N-1))/(12*dx)));
                        vynew(i,j) = vy(i,j) - dt*(vy(i,j)*(-vy(1,j)+8*vy(i+1,j)-8*vy(i-1,j)+vy(i-2,j))/(12*dy) + g*((-h(1,j)+8*h(i+1,j)-8*h(i-1,j)+h(i-2,j))/(12*dy)));
                        if h(i,j) > B(i,j)
                            hnew(i,j) = h(i,j) - dt*(vx(i,j)*((h(i,j+1)-h(i,N))/(2*dx)-dB_dx(i,j)) + vy(i,j)*((h(i+1,j)-h(i-1,j))/(2*dx)-dB_dy(i,j)));
                        end
                    elseif j == 2
                        vxnew(i,j) = vx(i,j) - dt*(vx(i,j)*(-vx(i,j+2)+8*vx(i,j+1)-8*vx(i,j-1)+vx(i,N))/(12*dx) + g*((-h(i,j+2)+8*h(i,j+1)-8*h(i,j-1)+h(i,N))/(12*dx)));
                        vynew(i,j) = vy(i,j) - dt*(vy(i,j)*(-vy(1,j)+8*vy(i+1,j)-8*vy(i-1,j)+vy(i-2,j))/(12*dy) + g*((-h(1,j)+8*h(i+1,j)-8*h(i-1,j)+h(i-2,j))/(12*dy)));
                        if h(i,j) > B(i,j)
                            hnew(i,j) = h(i,j) - dt*(vx(i,j)*((h(i,j+1)-h(i,j-1))/(2*dx)-dB_dx(i,j)) + vy(i,j)*((h(i+1,j)-h(i-1,j))/(2*dy) - dB_dy(i,j)));
                        end
                    elseif j == N-1
                        vxnew(i,j) = vx(i,j) - dt*(vx(i,j)*(-vx(i,1)+8*vx(i,j+1)-8*vx(i,j-1)+vx(i,j-2))/(12*dx) + g*((-h(i,1)+8*h(i,j+1)-8*h(i,j-1)+h(i,j-2))/(12*dx)));
                        vynew(i,j) = vy(i,j) - dt*(vy(i,j)*(-vy(1,j)+8*vy(i+1,j)-8*vy(i-1,j)+vy(i-2,j))/(12*dy) + g*((-h(1,j)+8*h(i+1,j)-8*h(i-1,j)+h(i-2,j))/(12*dy)));
                        if h(i,j) > B(i,j)
                            hnew(i,j) = h(i,j) - dt*(vx(i,j)*((h(i,j+1)-h(i,j-1))/(2*dx)-dB_dx(i,j)) + vy(i,j)*((h(i+1,j)-h(i-1,j))/(2*dy) - dB_dy(i,j)));
                        end
                    elseif j == N
                        vxnew(i,j) = vx(i,j) - dt*(vx(i,j)*(-vx(i,2)+8*vx(i,1)-8*vx(i,j-1)+vx(i,j-2))/(12*dx) + g*((-h(i,2)+8*h(i,1)-8*h(i,j-1)+h(i,j-2))/(12*dx)));
                        vynew(i,j) = vy(i,j) - dt*(vy(i,j)*(-vy(1,j)+8*vy(i+1,j)-8*vy(i-1,j)+vy(i-2,j))/(12*dy) + g*((-h(1,j)+8*h(i+1,j)-8*h(i-1,j)+h(i-2,j))/(12*dy)));
                        if h(i,j) > B(i,j)
                            hnew(i,j) = h(i,j) - dt*(vx(i,j)*((h(i,1)-h(i,j-1))/(2*dx)-dB_dx(i,j)) + vy(i,j)*((h(i+1,j)-h(i-1,j))/(2*dy) - dB_dy(i,j)));
                        end
                    end
                elseif i == N
                    if j == 1
                        vxnew(i,j) = vx(i,j) - dt*(vx(i,j)*(-vx(i,j+2)+8*vx(i,j+1)-8*vx(i,N)+vx(i,N-1))/(12*dx) + g*((-h(i,j+2)+8*h(i,j+1)-8*h(i,N)+h(i,N-1))/(12*dx)));
                        vynew(i,j) = vy(i,j) - dt*(vy(i,j)*(-vy(2,j)+8*vy(1,j)-8*vy(i-1,j)+vy(i-2,j))/(12*dy) + g*((-h(2,j)+8*h(1,j)-8*h(i-1,j)+h(i-2,j))/(12*dy)));
                        if h(i,j) > B(i,j)
                            hnew(i,j) = h(i,j) - dt*(vx(i,j)*((h(i,j+1)-h(i,N))/(2*dx)-dB_dx(i,j)) + vy(i,j)*((h(1,j)-h(i-1,j))/(2*dy) - dB_dy(i,j)));
                        end
                    elseif j == 2
                        vxnew(i,j) = vx(i,j) - dt*(vx(i,j)*(-vx(i,j+2)+8*vx(i,j+1)-8*vx(i,j-1)+vx(i,N))/(12*dx) + g*((-h(i,j+2)+8*h(i,j+1)-8*h(i,j-1)+h(i,N))/(12*dx)));
                        vynew(i,j) = vy(i,j) - dt*(vy(i,j)*(-vy(2,j)+8*vy(1,j)-8*vy(i-1,j)+vy(i-2,j))/(12*dy) + g*((-h(2,j)+8*h(1,j)-8*h(i-1,j)+h(i-2,j))/(12*dy)));
                        if h(i,j) > B(i,j)
                            hnew(i,j) = h(i,j) - dt*(vx(i,j)*((h(i,j+1)-h(i,j-1))/(2*dx)-dB_dx(i,j)) + vy(i,j)*((h(1,j)-h(i-1,j))/(2*dy) - dB_dy(i,j)));
                        end
                    elseif j == N-1
                        vxnew(i,j) = vx(i,j) - dt*(vx(i,j)*(-vx(i,1)+8*vx(i,j+1)-8*vx(i,j-1)+vx(i,j-2))/(12*dx) + g*((-h(i,1)+8*h(i,j+1)-8*h(i,j-1)+h(i,j-2))/(12*dx)));
                        vynew(i,j) = vy(i,j) - dt*(vy(i,j)*(-vy(2,j)+8*vy(1,j)-8*vy(i-1,j)+vy(i-2,j))/(12*dy) + g*((-h(2,j)+8*h(1,j)-8*h(i-1,j)+h(i-2,j))/(12*dy)));
                        if h(i,j) > B(i,j)
                            hnew(i,j) = h(i,j) - dt*(vx(i,j)*((h(i,j+1)-h(i,j-1))/(2*dx)-dB_dx(i,j)) + vy(i,j)*((h(1,j)-h(i-1,j))/(2*dy) - dB_dy(i,j)));
                        end
                    elseif j == N
                        vxnew(i,j) = vx(i,j) - dt*(vx(i,j)*(-vx(i,2)+8*vx(i,1)-8*vx(i,j-1)+vx(i,j-2))/(12*dx) + g*((-h(i,2)+8*h(i,1)-8*h(i,j-1)+h(i,j-2))/(12*dx)));
                        vynew(i,j) = vy(i,j) - dt*(vy(i,j)*(-vy(2,j)+8*vy(1,j)-8*vy(i-1,j)+vy(i-2,j))/(12*dy) + g*((-h(2,j)+8*h(1,j)-8*h(i-1,j)+h(i-2,j))/(12*dy)));
                        if h(i,j) > B(i,j)
                            hnew(i,j) = h(i,j) - dt*(vx(i,j)*((h(i,1)-h(i,j-1))/(2*dx)-dB_dx(i,j)) + vy(i,j)*((h(1,j)-h(i-1,j))/(2*dy) - dB_dy(i,j)));
                        end
                    end
                elseif j == 1
                        vxnew(i,j) = vx(i,j) - dt*(vx(i,j)*(-vx(i,j+2)+8*vx(i,j+1)-8*vx(i,N)+vx(i,N-1))/(12*dx) + g*((-h(i,j+2)+8*h(i,j+1)-8*h(i,N)+h(i,N-1))/(12*dx)));
                        vynew(i,j) = vy(i,j) - dt*(vy(i,j)*(-vy(i+2,j)+8*vy(i+1,j)-8*vy(i-1,j)+vy(i-2,j))/(12*dy) + g*((-h(i+2,j)+8*h(i+1,j)-8*h(i-1,j)+h(i-2,j))/(12*dy)));
                        if h(i,j) > B(i,j)
                            hnew(i,j) = h(i,j) - dt*(vx(i,j)*((h(i,j+1)-h(i,N))/(2*dx)-dB_dx(i,j)) + vy(i,j)*((h(i+1,j)-h(i-1,j))/(2*dy) - dB_dy(i,j)));
                        end
                elseif j == 2
                        vxnew(i,j) = vx(i,j) - dt*(vx(i,j)*(-vx(i,j+2)+8*vx(i,j+1)-8*vx(i,j-1)+vx(i,N))/(12*dx) + g*((-h(i,j+2)+8*h(i,j+1)-8*h(i,j-1)+h(i,N))/(12*dx)));
                        vynew(i,j) = vy(i,j) - dt*(vy(i,j)*(-vy(i+2,j)+8*vy(i+1,j)-8*vy(i-1,j)+vy(i-2,j))/(12*dy) + g*((-h(i+2,j)+8*h(i+1,j)-8*h(i-1,j)+h(i-2,j))/(12*dy)));
                        if h(i,j) > B(i,j)
                            hnew(i,j) = h(i,j) - dt*(vx(i,j)*((h(i,j+1)-h(i,j-1))/(2*dx)-dB_dx(i,j)) + vy(i,j)*((h(i+1,j)-h(i-1,j))/(2*dy) - dB_dy(i,j)));
                        end
                elseif j == N-1
                    vxnew(i,j) = vx(i,j) - dt*(vx(i,j)*(-vx(i,1)+8*vx(i,j+1)-8*vx(i,j-1)+vx(i,j-2))/(12*dx) + g*((-h(i,1)+8*h(i,j+1)-8*h(i,j-1)+h(i,j-2))/(12*dx)));
                    vynew(i,j) = vy(i,j) - dt*(vy(i,j)*(-vy(i+2,j)+8*vy(i+1,j)-8*vy(i-1,j)+vy(i-2,j))/(12*dy) + g*((-h(i+2,j)+8*h(i+1,j)-8*h(i-1,j)+h(i-2,j))/(12*dy)));
                    if h(i,j) > B(i,j)
                        hnew(i,j) = h(i,j) - dt*(vx(i,j)*((h(i,j+1)-h(i,j-1))/(2*dx)-dB_dx(i,j)) + vy(i,j)*((h(i+1,j)-h(i-1,j))/(2*dy)-dB_dy(i,j)));
                    end
                elseif j == N
                    vxnew(i,j) = vx(i,j) - dt*(vx(i,j)*(-vx(i,2)+8*vx(i,1)-8*vx(i,j-1)+vx(i,j-2))/(12*dx) + g*((-h(i,2)+8*h(i,1)-8*h(i,j-1)+h(i,j-2))/(12*dx)));
                    vynew(i,j) = vy(i,j) - dt*(vy(i,j)*(-vy(i+2,j)+8*vy(i+1,j)-8*vy(i-1,j)+vy(i-2,j))/(12*dy) + g*((-h(i+2,j)+8*h(i+1,j)-8*h(i-1,j)+h(i-2,j))/(12*dy)));
                    if h(i,j) > B(i,j)
                        hnew(i,j) = h(i,j) - dt*(vx(i,j)*((h(i,1)-h(i,j-1))/(2*dx)-dB_dx(i,j)) + vy(i,j)*((h(i+1,j)-h(i-1,j))/(2*dy) - dB_dy(i,j)));
                    end
                else %% se nao estiver a iterar na fronteira ou em R <= 10
                    vxnew(i,j) = vx(i,j) - dt*(vx(i,j)*(-vx(i,j+2)+8*vx(i,j+1)-8*vx(i,j-1)+vx(i,j-2))/(12*dx) + g*((-h(i,j+2)+8*h(i,j+1)-8*h(i,j-1)+h(i,j-2))/(12*dx)));
                    vynew(i,j) = vy(i,j) - dt*(vy(i,j)*(-vy(i+2,j)+8*vy(i+1,j)-8*vy(i-1,j)+vy(i-2,j))/(12*dy) + g*((-h(i+2,j)+8*h(i+1,j)-8*h(i-1,j)+h(i-2,j))/(12*dy)));
                    if h(i,j) > B(i,j)
                        hnew(i,j) = h(i,j) - dt*(vx(i,j)*((h(i,j+1)-h(i,j-1))/(2*dx)-dB_dx(i,j)) + vy(i,j)*((h(i+1,j)-h(i-1,j))/(2*dy)-dB_dy(i,j)));
                    end
                end
            end
        end
    end
    hmax(t_count) = max(max(h));
    t_hmax(t_count) = ti;
    vol(t_count) = volume(h,x,y);
    t_count = t_count+1;
    vx = vxnew;
    vy = vynew;
    h = hnew;
    str = sprintf('Tempo %f (s)', ti);
    figure(03);
    h2=surf(X, Y, h); %grafico da onda em x e y
    set(h2,'LineStyle','none');
    title(str);
    xlabel('x(m)');
    ylabel('y(m)');
    zlabel('Wave Height (m)');
    ti = ti + dt;
    pause(0.01)
end
figure(04);
plot(t_hmax,hmax);title ('Altura maxima da onda em funcao do tempo');
figure(05);
plot(t_hmax,vol); title('Estimativa do volume de agua deslocado');
            
                        
                        
                    
                    