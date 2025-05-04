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
    .container{
        display: flex;
    }
    .col{
        flex: 1;
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



## Use the Jacobi Method
<div class="container">
<div class="col">

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
</div>


<div class="col">

Exact solution
$$
x^*
$$
Error
$$e_0 = x^*-x_0$$
Residual
$$r_0 = b-Ax_0 = A e_0$$

</div>
</div>

---

## For Model Problem

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


As matrix

$$
A=\frac{1}{h^2}\left[\begin{array}{ccc}
2 & -1 & \\
-1 & \ddots & \ddots \\
& \ddots &
\end{array}\right]
$$

---

<!-- _header: Model Problem -->
<div class="container">
<div class="col">

Look at the matrix
$$
A=\left[\begin{array}{ccccc}
2 & -1 & & & \\
-1 & 2 & -1 & & \\
& \ddots & \ddots & \ddots & \\
& & -1 & 2 & -1 \\
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

</div>
<div class="col">

Eigenvectors for $n=64$

<style scoped>
img[alt~="eigen"]{
  display: block;
  margin: 0 auto;
  width: 70%;
}
</style>
![eigen](./fig/model_eigen.svg)

Performs like $\textcolor{#2ba02b}{\text{low}}$ and $\textcolor{#ff7f0f}{\text{high}}$ frequencies

</div>
</div>

---

## What does Jacobi do to error?

<div class="container">

<div class="col">

The error propagation
$$
e \leftarrow T e \quad T=I-D^{-1} A
$$

$$
T=\left[\begin{array}{ccccc}
0 & 1/2 & & & \\
1/2 & 0 & 1/2 & & \\
& \ddots & \ddots & \ddots & \\
& & 1/2 & 0 & 1/2 \\
& & & 1/2 & 0
\end{array}\right]
$$

It is averaging

$$
\textcolor{#2ba02b}{e_i^{\text {new }}} \leftarrow \frac{1}{2}\left(\textcolor{#ff7f0f}{e_{i-1}^{\text {old }}}+\textcolor{#ff7f0f}{e_{i+1}^{\text {old }}}\right)
$$

</div>

<div class="col">

For different types of error

<style scoped>
img[alt~="avg"]{
  display: block;
  margin: 0 auto;
  width: 60%;
}
</style>
![avg](./fig/avg.svg)

It *averages out* certain frequency quickly

</div>

</div>


---

## From Jacobi to weighted-Jacobi

<div class="container">

<style scoped>
img[alt~="eigen"]{
  display: block;
  margin: 0 auto;
  width: 90%;
}
</style>

<div class="col">

$$
u \leftarrow u+D^{-1} r
$$

![eigen](./fig/wj1.svg)

</div>

<div class="col">

$$
u \leftarrow u+\omega D^{-1} r , \omega = 2 / 3
$$

![eigen](./fig/wj2.svg)

</div>

</div>

Better but why $2/3$?

---

<!-- _header: Fourier Analysis of Errors -->

Using the eigenvectors of the error propagation $T$ as a basis for the error space
$$
e_0=\sum_{k=1}^n c_k v_k
$$

Then the error transforms like
$$
e \leftarrow\left(I-\omega D^{-1} A\right) e=T e
$$

$$
\begin{aligned}
e_1 = T e_0 & =\sum_{k=1}^n c_k T v_k \\
& =\sum_{k=1}^n c_k \lambda_k v_k
\end{aligned}
$$

Error on the direction of the $v_k$, or frequency $k$, is reduced by the magnitude of $\lambda_k$

---

<!-- _header: Fourier Analysis of Errors -->

Eigenvalues of $T$ corresponding to different $\omega$

<style scoped>
img[alt~="eigen"]{
  display: block;
  margin: 0 auto;
  width: 55%;
}
</style>
![eigen](./fig/omega_eigen.svg)

Weighted Jacobi with $\omega = 2/3$ dampens errors in high-frequency modes

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
