# Projeto - Pacote de onda Quântico

# Bibliotecas
import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import fft, ifft
from matplotlib.animation import FuncAnimation

# Constantes em um novo sistema de medidas
hbar = 1.0  # Constante de Planck reduzida
m = 1.0     # Massa da partícula

# Grade espacial e temporal (array 1D)
N = 1000
L = 100.0 # tamanho da grade
x = np.linspace(-L/2, L/2, N)
dx = x[1] - x[0] # distância entre os pontos
dt = 0.01 

# Criando a função de onda inicial psi(x,0)\pacote gaussiano
x0 = -L/4 #posição inicial
largura = L/20 #largura do pacote
k0 = 10.0 #momento inicial para controlar a velocidade 

#construção da função de onda
parte_gaussiana = np.exp(-(x-x0)**2 / (2*largura**2))
parte_onda_plana = np.exp(1j*k0*x) # 1j representa a parte imaginária
psi = parte_gaussiana*parte_onda_plana
psi = psi/np.sqrt(np.sum(np.abs(psi)**2)*dx) #normalização

#Gráfico densidade de probabilidade
'''densidade_prob = np.abs(psi)**2
plt. figure(figsize=(10,6))
plt.plot(x, densidade_prob, label = r'$|\psi(x,0)|^2$')
plt.title('Densidade de Probabilidade Inicial do Pacote de Onda')
plt.xlabel('Posição (x)')
plt.ylabel('Probabilidade')
plt.grid(True)
plt.legend()
plt.show()'''


# Resolvendo a equação de Schrodinguer no tempo - split step fourier
#dividimos a solução da equação em: cinética e potencial

# Criando o espaço dos momentos (k)
dk = 2*np.pi/L
k = np.hstack([np.arange(0,N/2), np.arange(-N/2, 0)])*dk

#Função da evolução temporal
def evolucao_passo_unico(psi, V, k, hbar, m, dt):
    #1 evolução por meio passo com potencial V
    psi = np.exp(-1j*V*dt/(2*hbar))*psi
    #2 mudança para o espaço dos momentos
    psi_k = fft(psi)
    psi_k = np.exp(-1j*hbar*k**2*dt/(2*m))*psi_k
    #3 voltando para o espaço da posição (FFT inversa)
    psi = ifft(psi_k)
    #4 evoluir pela segunda metade do passo com potencial V
    psi = np.exp(-1j*V*dt/ (2*hbar))*psi
    return psi

# Construindo uma barreira potencial
V0 = 50.0 #altura da barreira de potencial
largura_barreira = L / 50
V = np.zeros_like(x)
indices_barreira = np.abs(x) < (largura_barreira / 2) #define a região da barreira
V[indices_barreira] = V0

#---- Plotando estado inicial e a barreira ----
densidade_prob_inicial = np.abs(psi)**2

fig, ax1 = plt.subplots(figsize = (10, 6))
ax1.plot(x, densidade_prob_inicial, 'b-', label = r'$|\psi(x,0)|^2$')
ax1.set_xlabel('Posição (x)')
ax1.set_ylabel('Probabilidade', color = 'b')
ax1.tick_params('y', color = 'b')
ax1.grid(True)

ax2 = ax1.twinx()
ax2.plot(x, V, 'r--', label = 'Potencial V(x)')
ax2.set_ylabel('Energia', color = 'r')
ax2.tick_params('y', colors = 'r')

plt.title('Estado Inicial e Barreira de Potencial')
fig.tight_layout()
#plt.show()
# ----------------------------------------

# Loop da simulação
numero_passos = 800
passos_plot = 10

psi_resultados = [] # guarda os quadros da animação
psi_atual = psi.copy() # começa com a função de onda inicial

for i in range(numero_passos):
    psi_atual = evolucao_passo_unico(psi_atual, V, k, hbar, m, dt)
    #guarda o resultado a cada passo_plot
    if i % passos_plot ==0:
        psi_resultados.append(psi_atual.copy())
print(f'Simulação concluida. {len(psi_resultados)} quadros gravados para a animação.')

#------ Animação ------
fig, ax = plt.subplots(figsize = (10, 6))

linha_prob, = ax.plot([], [], 'b-', lw = 2)
ax.plot(x, V, 'r--', label = 'Potencial V(x)')
ax.set_xlim(-L/2, L/2)
ax.set_ylim(0, np.max(np.abs(psi)**2)* 1.1)
ax.set_xlabel('Posição(x)')
ax.set_ylabel('Probabilidade')
ax.legend()
ax.grid(True)

def update(frame):
    psi_t = psi_resultados[frame]
    prob_densidade_t = np.abs(psi_t)**2
    linha_prob.set_data(x, prob_densidade_t)
    return linha_prob,

anim = FuncAnimation(fig, update, frames = len (psi_resultados), blit = True, interval = 20)
plt.show()