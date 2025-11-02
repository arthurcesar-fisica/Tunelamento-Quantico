#  Simulação de Tunelamento Quântico com Método Split-Step Fourier

Este projeto é uma simulação em Python da dinâmica de um pacote de onda quântico unidimensional, resolvendo a **Equação de Schrödinger Dependente do Tempo (TDSE)**.

O objetivo é visualizar o fenômeno do **tunelamento quântico**: um pacote de onda (representando uma partícula) incide sobre uma barreira de potencial finita. Mesmo com energia cinética *menor* que a altura da barreira, uma parte da função de onda "atravessa" a barreira, resultando em uma probabilidade não nula de encontrar a partícula do outro lado.

A simulação utiliza o eficiente **Método Split-Step Fourier** para evoluir a função de onda no tempo e `matplotlib` para animar a densidade de probabilidade $|\psi(x,t)|^2$.

## Animação do Resultado

O script gera uma animação que mostra o pacote de onda Gaussiano se aproximando, colidindo e se dividindo ao interagir com a barreira de potencial.


`![Animação do Tunelamento](./tunelamento_anim.gif)`

---

###  Contexto Físico: A Equação de Schrödinger

O nosso objetivo é resolver a TDSE 1D:
$$i\hbar \frac{\partial \psi(x,t)}{\partial t} = \hat{H} \psi(x,t) = \left( \hat{T} + \hat{V} \right) \psi(x,t)$$

Onde $\hat{T} = -\frac{\hbar^2}{2m}\frac{\partial^2}{\partial x^2}$ é o operador de energia cinética e $\hat{V} = V(x)$ é o operador de energia potencial.

A solução formal seria $\psi(t+\Delta t) = e^{-i\hat{H}\Delta t / \hbar} \psi(t)$. O problema é que os operadores $\hat{T}$ e $\hat{V}$ não comutam, então não podemos simplesmente separar $e^{-i(\hat{T}+\hat{V})\Delta t / \hbar}$ em $e^{-i\hat{T}\Delta t / \hbar} \cdot e^{-i\hat{V}\Delta t / \hbar}$.

---

### 💻 A Técnica: O Método Split-Step Fourier

Aqui entra o "pulo do gato" computacional. Usamos a aproximação de Trotter (ou Baker-Campbell-Hausdorff) para "dividir" (split) a evolução em pequenos passos:

$$\psi(t+\Delta t) \approx e^{-i\hat{V} \Delta t / 2\hbar} \cdot e^{-i\hat{T} \Delta t / \hbar} \cdot e^{-i\hat{V} \Delta t / 2\hbar} \psi(t)$$

Esta aproximação é muito precisa para $\Delta t$ pequenos. A sua grande vantagem é:

1.  **Parte do Potencial ($e^{-i\hat{V} \Delta t / 2\hbar}$):** Este operador é diagonal no **espaço de posição (x)**. A evolução é uma simples multiplicação: `psi = np.exp(-1j * V * dt / (2*hbar)) * psi`.
2.  **Parte Cinética ($e^{-i\hat{T} \Delta t / \hbar}$):** Este operador é diagonal no **espaço de momento (k)**, onde $\hat{T} \rightarrow \frac{\hbar^2 k^2}{2m}$.

O **Método Fourier** é a ponte entre esses dois espaços. A **Transformada Rápida de Fourier (FFT)** nos permite pular de $x$ para $k$ (e voltar com a `ifft`) de forma extremamente rápida (complexidade $\mathcal{O}(N \log N)$).

#### O Algoritmo (implementado em `evolucao_passo_unico`):
1.  Evolui a função de onda por meio passo de potencial $\Delta t/2$ no espaço $x$.
2.  Usa `fft(psi)` para ir para o espaço $k$.
3.  Evolui um passo completo de cinética $\Delta t$ no espaço $k$.
4.  Usa `ifft(psi_k)` para voltar para o espaço $x$.
5.  Evolui a segunda metade do passo de potencial $\Delta t/2$ no espaço $x$.

---

## 🚀 Como Executar

O script foi escrito em Python e depende de algumas bibliotecas científicas padrão.

**1. Pré-requisitos**

Instale as bibliotecas necessárias:
```bash
pip install numpy matplotlib