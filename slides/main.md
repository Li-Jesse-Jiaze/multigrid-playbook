---
marp: true
paginate: true
math: katex
style: |
  section .sidebox{
    position: absolute;
    bottom: 2rem;
    right: 1rem;
    width: 30%;
    padding: 0.5rem .5rem;
    border: 0px;
    border-radius: 6px;
    background: rgba(0,0,0,.1);
    font-size: 0.85em;
  }
    .container{
        display: flex;
    }
    .col{
        flex: 1;
    }
    .site {
    color: gray;
    font-size: 0.5em;
    margin-top: -1em;
    display: block;
    }
---

<style scoped>
img[alt~="cover"]{
  position:absolute; right:1rem; bottom:2rem; width:70%;
}
</style>
![cover](./fig/cover.svg)

# Multigrid Methods

Y. Saad, *Iterative Methods for Sparse Linear Systems*, §§13.3–13.5

**Presenter: Jiaze Li**

<br>
<br>
<br>
<br>

---

## Model Problem

1D model problem

$$
\begin{aligned}
-u_{x x} & =f \\
u(0) & =u(1)=0
\end{aligned}
$$


With finite differences

$$
\frac{-u_{i-1}+2 u_i-u_{i+1}}{h^2}=f_i \quad i=1, \ldots, n \quad u_0=u_{n+1}=0
$$


As matrix form

$$
\frac{1}{h^2}\left[\begin{array}{ccc}
2 & -1 & \\
-1 & \ddots & \ddots \\
& \ddots &
\end{array}\right] \left[\begin{array}{c}
u_1 \\
\vdots \\
u_{n}
\end{array}\right] = 
\left[\begin{array}{c}
f_1 \\
\vdots \\
f_{n}
\end{array}\right]
$$

---

### Use the Jacobi Method
<br>

Solve the problem
$$
A x=b
$$

Iteratively
$$
x^{(1)}=x^{(0)}+D^{-1} r^{(0)}
$$

<br>


> Or think of it as a gradient descent with
> fixed step size $1$ and preconditioner $D^{-1}$

<div class="sidebox">

Exact solution
$$
x^*
$$
Error
$$e^{(0)} = x^*-x^{(0)}$$
Residual
$$r^{(0)} = b-Ax^{(0)} = A e^{(0)}$$

</div>

---

## Some Convergence Test

The number of iterations required to get the tolerance of 1.0E-10 (RMS residual norm)

<style scoped>
img[alt~="n_iter"]{
  display: block;
  margin: 0 auto;
  width: 40%;
}
</style>
![n_iter](./fig/n_iter.svg)

<span class="site">Burkardt, J., 2011. Jacobi Iterative Solution of Poisson’s Equation in 1D [online]</span>

The computational cost is about $\mathcal{O}(n^3)$. But why is that? Let's look at the details.

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

The eigenvectors (are Fourier modes)

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

### What does Jacobi do to error?

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

It is averaging (mean filter)

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

### From Jacobi to weighted-Jacobi

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

Like a low-pass filter. Why?

---

### Fourier Analysis

Using the eigenvectors of the error propagation $T$ as a basis for the error space
$$
e^{(0)}=\sum_{k=1}^n c_k v_k
$$

Then the error transforms like
$$
e \leftarrow\left(I-\omega D^{-1} A\right) e=T e
$$

$$
e^{(\sigma)} = T e^{(0)} =\sum_{k=1}^n c_k T^\sigma v_k =\sum_{k=1}^n c_k \lambda_k^\sigma v_k
$$

Converge with the factor $\|T\|$ or $\rho(T)$

---

<!-- _header: Fourier Analysis -->

For Jacobi and weighted-Jacobi
$$
\begin{aligned}
T=I-D^{-1} A & \longrightarrow \lambda_k=1-\frac{1}{2} \cdot 4 \cdot \sin ^2\left(\frac{k \pi}{2(n+1)}\right) \\
T=I-(2 / 3) D^{-1} A & \longrightarrow \lambda_k=1-\frac{2}{3} \cdot \frac{1}{2} \cdot 4 \cdot \sin ^2\left(\frac{k \pi}{2(n+1)}\right)
\end{aligned}
$$

The spectral radius $\rho_n \approx 1-{C}/{n^2}$ (confirms the quadratic **#iteration**)

But
$$
e^{(\sigma)} =\sum_{k=1}^n c_k \lambda_k^\sigma v_k
$$
error on different direction $v_k$, or frequency $k$, is reduced by different magnitude of $\lambda_k$


---

<!-- _header: Fourier Analysis -->

$\lambda_k$ corresponding to different $\omega$

<style scoped>
img[alt~="eigen"]{
  display: block;
  margin: 0 auto;
  width: 40%;
}
</style>
![eigen](./fig/omega_eigen.svg)

For $\omega = 2/3$, the high frequency part of the error is reduced by (at least) $1/3$ 

> **Smoothing factor** of $T$: the maximum magnitude of the upper half spectrum

Maybe it's possible to apply this 1/3 reduction on half the error, on and on?

---

### How can we turn low frequency errors to high frequency?

Sampling on a coarse grid

<style scoped>
img[alt~="eigen"]{
  display: block;
  margin: 0 auto;
  width: 75%;
}
</style>
![eigen](./fig/coarse_grid.svg)

After this, smooth modes look like oscillatory modes, with little loss of information.

---

### So far we got

Classical iterative methods converge slowly as $n \uparrow$, but have the smoothing property.

Low-frequency information is well approximated on the coarse grid and becomes high-frequency.

### Next questions

How to transfer between fine and coarse?

What do we do “solve” on a coarse grid?

---


### Remind the Projection Methods

Look for the "best" update:
$$
x^{(1)} \leftarrow x^{(0)}+u
$$

Over a smaller space (with $A$-norm)
$$
\min _{u \in \text { span }\{V\}}\left\|x^*-x^{(1)}\right\|_A
$$

So the update looks like
$$
x^{(1)}=x^{(0)}+V\left(V^\top A V\right)^{-1} V^\top r^{(0)}
$$

The error
$$
e^{(1)} = (I-\underbrace{\textcolor{green}{V\left(V^\top A V\right)^{-1} V^\top A}}_\text{A-orthogonal projection}) e^{(0)}
$$

---

### Prolongation: From Coarse to Fine

<style scoped>
img[alt~="c2f"]{
  display: block;
  margin: 0 auto;
  width: 60%;
}
</style>
![c2f](./fig/c2f.SVG)

Construct an operator
$$
P: \Omega^{2 h} \rightarrow \Omega^h
$$
where $\Omega^{2 h}$ is the coarse grid and $\Omega^{h}$ is the fine grid

---

### Interpolation

<div class="container">

<div class="col">

$$
\begin{aligned}
v_{2 i}^h & =v_i^{2 h}\\
v_{2 i+1}^h & =\frac{1}{2}\left(v_i^{2 h}+v_{i+1}^{2 h}\right)
\end{aligned}
$$

Or in matrix form
$$
P = \frac{1}{2}\left[\begin{array}{cccc}
1 & & & \\
2 & & & \\
1 & 1 & & \\
& 2 & & \\
& & \vdots & \\
& & 1 & 1 \\
& & & 2 \\
& & & 1
\end{array}\right]
$$

</div>

<div class="col">

Looks like
<style scoped>
img[alt~="interp"]{
  display: block;
  margin: 0 auto;
  width: 80%;
}
</style>
![interp](./fig/interp.svg)

Notice: $P$ is full rank

</div>

</div>

---

In projection method
$$
x^{(1)}=x^{(0)}+V\left(V^\top A V\right)^{-1} V^\top r^{(0)}
$$

If $V$ is $\text { span }\{p\}$, or just $P$

$$
x^{(1)}=x^{(0)}+P\left(P^\top A P\right)^{-1} P^\top r^{(0)}
$$

We get the **two-grid method**

---

### Two-grid Method

$$
x^{(1)}=x^{(0)}+P\left(P^\top A P\right)^{-1} P^\top r^{(0)}
$$

$$
\begin{array}{lc}
\text {1. Given } & x^{(0)} \\
\text {2. Smooth a few times } & x^{(0)} \leftarrow x^{(0)}+\omega D^{-1} A r^{(0)} \\
\text {3. Form residual } & r^{(0)}=b-A x^{(0)} \\
\text {4. Restrict the residual } & P^\top r^{(0)} \\
\text {5. Solve the coarse problem } & P^\top A P \delta=P^\top r^{(0)} \quad \\
\text {6. Interpolate the approx error } & P \delta \\
\text {7. Correct} & x^{(1)}= x^{(0)}+P \delta
\end{array}
$$

$R = P^\top$ is the **restriction** here
$A_c = P^\top A P$ is the **coarse level operator** here

---

Let's look at the coarse level operator $A_c$

<style scoped>
img[alt~="glk"]{
  display: block;
  margin: 0 auto;
  width: 80%;
}
</style>
![glk](./fig/galerkin.svg)

<div class="container">

<div class="col">


The restriction $R$ don't have to be $P^\top$

See if

$$
R = P^\top = \frac{1}{2}\left[\begin{array}{cccccccccc}
1 & 2 & 1 & & & & & & & \\
& & 1 & 2 & 1 & & & & & \\
& & & & \cdots & \cdots & \cdots & & & \\
& & & & & & & 1 & 2 & 1
\end{array}\right]
$$

</div>

<div class="col">

<style scoped>
img[alt~="restrict"]{
  display: block;
  margin: 0 auto;
  width: 80%;
}
</style>
![restrict](./fig/restrict.svg)

</div>

</div>

---

### Other options for restriction

Make sure that each row of $R$ sums to 1

<div class="container">

<style scoped>
img[alt~="restrict"]{
  display: block;
  margin: 0 auto;
  width: 80%;
}
</style>

<div class="col">


Weighted $P^\top$

$$
\frac{1}{4}\left[\begin{array}{cccccccccc}
1 & 2 & 1 & & & & & & & \\
& & 1 & 2 & 1 & & & & & \\
& & & & \cdots & \cdots & \cdots & & & \\
& & & & & & & 1 & 2 & 1
\end{array}\right]
$$

![restrict](./fig/restrict_weighted.svg)

</div>

$\quad$

<div class="col">

Injection

$$
\left[\begin{array}{cccccccccc}
0 & 1 & 0 & & & & & & & \\
& & 0 & 1 & 0 & & & & & \\
& & & & \cdots & \cdots & \cdots & & & \\
& & & & & & & 0 & 1 & 0
\end{array}\right]
$$

<style scoped>
img[alt~="restrict"]{
  display: block;
  margin: 0 auto;
  width: 80%;
}
</style>
![restrict](./fig/restrict_injection.svg)

</div>

</div>

---

### Algorithm: Two-grid Method

<div class="container">

<div class="col">

1. Smooth $\nu_{p r e}$ times on $A u=f$
2. Compute $r=f-A u$
3. Compute $r_c=R r$
4. Solve $A_c \delta_c=r_c$
5. Interpolate $\hat{\delta}=P \delta_c$
6. Correct $u \leftarrow u+\hat{\delta}$
7. Smooth $\nu_{p o s t}$ times on $A u=f$

</div>

<div class="col">

<style scoped>
img[alt~="vcycle"]{
  display: block;
  margin: 0 auto;
  width: 80%;
}
</style>
![vcycle](./fig/vcycle.svg)

<center>A two-level V cycle</center>

</div>

</div>

---

### How Accurate is Multigrid?

$u^*-u_h^* \quad$ Discretization error
$u_h^*-u_h \quad$ Algebraic error

TODO: Show that the total error is limited by the discretization error

---

### Convergence of Two-grid Method

TODO: Analyze the convergence via Fourier and energy norm (by writing) ⭐

---

### Multigrid

Do two-grid recursively

TODO: W/V cycle
