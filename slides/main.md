---
marp: true
paginate: true
math: katex
style: |              # 全局样式（作用于所有幻灯片）
  /* 右侧信息框的通用样式 */
  section .sidebox{
    position: absolute;   /* 关键：绝对定位 */
    bottom: 2rem;            /* 距离顶部，可自行调整 */
    right: 1rem;        /* 与右边保持间距 */
    width: 30%;           /* 盒子宽度 */
    padding: 0.5rem .5rem;
    border: 0px;
    border-radius: 6px;
    background: rgba(0,0,0,.1);
    font-size: 0.85em;    /* 字号稍微缩小一点 */
  }
---

<style scoped>
img[alt~="cover"]{
  position:absolute; right:1rem; bottom:2rem; width:70%;
}
</style>
![cover](./fig/cover.svg)

# Multigrid Methods

> Y. Saad, *Iterative Methods for Sparse Linear Systems*, §§13.3–13.5

Presenter: Jiaze Li

<br>
<br>
<br>
<br>

---

## Start with the Jacobi Method

Solve the problem
$$
A x=b
$$

Iteratively
$$
x_1 \leftarrow x_0 + \text { update }
$$

Ideally
$$
x_1=x_0+e_0=x_0+A^{-1} r_0
$$

Not practical so
$$
x_1=x_0+D^{-1} r_0
$$


<div class="sidebox">

Exact solution: $x^*$</li>
Error: $e_0 = x^*-x_0$</li>
Residual: $r_0 = b-Ax_0 = A e_0$

</div>

---

## Model Problem

1D model problem

$$
\begin{aligned}
-u_{x x} & =f \\
u(0) & =u(1)=0
\end{aligned}
$$


Finite differences

$$
\frac{-u_{i-1}+2 u_i-u_{i+1}}{h^2}=f_i \quad i=1, \ldots, n \quad u_0=u_{n+1}=0
$$


A model matrix problem

$$
A=\frac{1}{h^2}\left[\begin{array}{ccc}
2 & -1 & \\
-1 & \ddots & \ddots \\
& \ddots &
\end{array}\right]
$$

---

<!-- _header: Model Problem -->


Look at the matrix
$$
A=\left[\begin{array}{rrrrr}
2 & -1 & & & \\
-1 & 2 & -1 & & \\
& -1 & 2 & -1 & \\
& & \ddots & & \\
& & & -1 & 2
\end{array}\right]
$$

The eigenvalues

$$
\lambda_k=4 \sin ^2\left(\frac{k \pi}{2(n+1)}\right)
$$

The eigenvectors

$$
\left(v_k\right)_j=\sin \left(\frac{(j+1) * k \pi}{n+1}\right)
$$

---

TODO:
- Show that error damped rapidly on some direction (Fourier)
- Why using 2/3 for weighted Jacobi
- Lead to MG Step #1: smoother
- Show the relation between Coarse-Fine and Smooth-Oscillatory
- Lead to transfer between fine and coarse

---

## Questions to resolve...

How to transfer between fine and coarse?

What do we do “solve” on a coarse grid?

---


## Remind the Projection Methods

Look for the "best" update:
$$
x_1 \leftarrow x_0+u
$$

Over a smaller space
$$
\min _{u \in \operatorname{span}\{V\}}\left\|x^*-x_1\right\|
$$

Then $u=V y$
$$
V^T V y=V^T e_0
$$

So the update looks like
$$
x_1=x_0+V\left(V^T V\right)^{-1} V^T e_0
$$

---

<!-- _header: Remind the Projection Methods -->

Look at the $A$-norm
$$
\min _{u \in \text { span }\{V\}}\left\|x^*-x_1\right\|_A
$$

Then
$$
V^T A V y=V^T A e_0 =V^T r_0
$$

So that
$$
x_1=x_0+V\left(V^T A V\right)^{-1} V^T r_0
$$

The error
$$
e_1 = (I-\underbrace{\textcolor{green}{V\left(V^T A V\right)^{-1} V^T A}}_\text{A-orthogonal projection}) e_0
$$

---

TODO:
- Restriction and Prolongation
- Two-grid
- How Accurate is it?
- Convergence
- Multigrid W/V cycle
