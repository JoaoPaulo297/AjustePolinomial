# Ajuste Polinomial por Mínimos Quadrados

Projeto em C para calcular o ajuste polinomial (regressão polinomial) por mínimos quadrados.  
O projeto monta a matriz de Vandermonde a partir de pontos \((x, y)\), resolve o sistema normal
\((A^T A) P = A^T Y\) e obtém os coeficientes `P` do polinômio de grau `k` que minimiza o erro quadrático.

---

## Funcionalidades

- Geração da matriz de Vandermonde a partir de dados.
- Cálculo de \(A^T A\) e \(A^T Y\).
- Inversão de matrizes por Gauss-Jordan (com pivotamento parcial).
- Cálculo dos coeficientes do polinômio: \(P = (A^T A)^{-1} A^T Y\).
- Cálculo do SSE (soma dos quadrados dos erros).
- Medição do tempo de execução do ajuste.
