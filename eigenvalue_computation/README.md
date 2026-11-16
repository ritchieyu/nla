# Eigenvalue Computation

References:

- CME 302 notes by Prof. Eric Darve
- Textbook by Trefethen and Bau

## Main Ideas

If $A\in\mathbb{R}^{n\times n}$, we want to find eigenvalues $\lambda$ and eigenvectors $v$ (with $v\neq 0$) such that $Av=\lambda v$. We approach this problem by considering a sequence of smaller eigenvalue problems. This style appears in other iterative approaches to solve other matrix problems involving $A$, as well.

### How much work?

Generally speaking,

```math
\textrm{work} = \textrm{steps} \times \textrm{work / step}
```

- For **direct methods**, there is often a sequence of $O(m)$ steps that, once executed, will give the desired result.

- For **iterative methods**, we care about **convergence**.
  The number of steps to converge (reach a desired accuracy) depends on spectral properties of $A$ (e.g., how eigenvalues are clustered).

In both cases, the work per step depends on matrix structure and how it is exploited by the algorithm.

### Other Remarks

- Direct vs. Iterative methods
  - In both cases, the best we can hope for is an **answer accurate to machine precision**
  - Both are capable of this

- Some direct methods beat the typical $\mathcal{O}(m^3)$ cost, but they often are not practical for two reasons:
  - They may lack stability
  - The hidden constant may be (galactically) huge

## The Krylov Subspace and other relevant objects

Let $A\in\mathbb{R}^{n\times n}$ and $q_1$ a unit vector.
The Krylov subspace $K_k(A,q_1)$ is

```math
K_k(A,q_1)\coloneqq\textrm{span}(q_1,Aq_1,\dots, A^{k-1}q_1)
```

If $K_k(A,q_1)$ has dimension $k$, then $q_1,Aq_1,\dots,A^{k-1}q_1$ must be linearly independent.
Further, there exists an orthonormal basis $q_1,\dots,q_k$ of $K_k(A,q_1)$.
Define

```math
Q_k\coloneqq\Big(q_1, q_2,\dots,q_k\Big)
```

Using $Q_k$, let $x\in K_k(A,q_1)$ and express $x=Q_ky$.
In general, $Ax$ may not remain in the Krylov subspace, but we can **project** it back:

```math
K_k(A,q_1)\ni (QQ_k^T)Ax=Q(Q_k^TAQ_k)x
```

We define

```math
H_k\coloneqq Q_k^TAQ_k
```

What does $H_k$ represent? Since $(QQ_k)^TAQ_ky=Q_kH_ky$, then $H_k$ represents the **projected action** of $A$, written in coordinates of $Q_k$.

## Arnoldi Iteration