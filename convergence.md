## Problem

Consider a matrix problem (s.p.d.) of the form
$$
A \boldsymbol{u}=\boldsymbol{f}, \quad A \in \mathbb{R}^{n \times n}
$$
Suppose we have a multilevel iteration process $\mathcal{M}$

$$
I-\mathcal{M} A=\left(I-M^\top A\right)^{\nu_{\text {pre }}}\left(I-P\left(P^\top A P\right)^{-1} P^\top A\right)(I-M A)^{\nu_{\text {post }}}
$$

so that

$$
\boldsymbol{e} \leftarrow(I-\mathcal{M} A) \boldsymbol{e}
$$
The iteration converges for any $\boldsymbol{f}$ and $\boldsymbol{u}_0$ iff $\rho(I-\mathcal{M} A)<1$.

we are seeking a method that yields **a bound on the error reduction in each iteration that is independent of** $n$.

## Notation

Fine grid $\Omega=\{1, \ldots, n\}=C \cup F$ and coarse grid $\Omega_c=C$. Interpolation / restriction:
$$
P: \Omega_c \rightarrow \Omega \quad \text { and } \quad R: \Omega \rightarrow \Omega_c
$$

$A$ is s.p.d., $D=\operatorname{diag}(A)-$ defining an inner product:

$$
\begin{aligned}
\langle\boldsymbol{u}, \boldsymbol{v}\rangle_A & =\langle A \boldsymbol{u}, \boldsymbol{v}\rangle \\
\langle\boldsymbol{u}, \boldsymbol{v}\rangle_D & =\langle D \boldsymbol{u}, \boldsymbol{v}\rangle \\
\langle\boldsymbol{u}, \boldsymbol{v}\rangle_{A D^{-1} A} & =\left\langle D^{-1} A \boldsymbol{u}, A \boldsymbol{v}\right\rangle
\end{aligned}
$$

Define $S$ as post-relaxation
(and $\mathcal{S}$ as the affine version $\boldsymbol{u} \leftarrow \mathcal{S}(\boldsymbol{u}, \boldsymbol{f})$ ):
$$
\boldsymbol{u} \leftarrow S \boldsymbol{u}+(I-S) A^{-1} \boldsymbol{f} \quad \text { or } \quad \boldsymbol{e} \leftarrow S \boldsymbol{e}
$$
## Algorithm

$$
\begin{aligned}
& \boldsymbol{u} \leftarrow \hat{\mathcal{S}}(\boldsymbol{u}, \boldsymbol{f}) \\
& \boldsymbol{r}_c \leftarrow R \boldsymbol{r} \\
& \boldsymbol{e}_c \leftarrow A_c^{-1} \boldsymbol{r}_c \\
& \hat{\boldsymbol{e}} \leftarrow P \boldsymbol{e}_c \\
& \boldsymbol{u} \leftarrow \boldsymbol{u}+\hat{\boldsymbol{e}} \\
& \boldsymbol{u} \leftarrow \mathcal{S}(\boldsymbol{u}, \boldsymbol{f})
\end{aligned}
$$
In general, consider relaxation as

$$
S=I-M A
$$
Assumptions:

- $M$ is norm convergent (in $A$ ): $\|S\|_A<1$
- $P$ is full rank
- $A$ is s.p.d.

Let the coarse grid correction step be

$$
T=I-P\left(P^\top A P\right)^{-1} P^\top A
$$
$T$ is an $A$-orthogonal projection onto the range of $P$: After coarse grid correction, the error is minimized in the energy norm over $\text{range}(P)$

## Focus on $V(0, 1)$

The $A$-adjoint of $S T$ is

$$
T S^{+} \quad S^{+}=I-M^\top A
$$


The symmetric $V(1,1)$ cycle is

$$
\begin{aligned}
(I-M A)\left(I-P\left(P^T A P\right)^{-1} P^T A\right)\left(I-M^T A\right) & =S T S^{+} \\
& =S T T S^{+}
\end{aligned}
$$


Since $\|S T\|_A=\left\|T S^{+}\right\|_A$ ($A$-adjoints) we have

$$
\left\|S T S^{+}\right\|_A=\|S T\|_A^2
$$


Ok, so we can focus focus on the $V(0,1)$ cycle, the other cycles follow.

## Upper Bound

We will measure convergence or reduction in $\|\boldsymbol{e}\|_A$. Note:

$$
\|\boldsymbol{e}\|_A^2=\|(I-T) \boldsymbol{e}\|_A^2+\|T \boldsymbol{e}\|_A^2
$$


For a $V(0,1)$ cycle, the reduction in $\boldsymbol{e}$ is

$$
\|S T e\|_A^2 \leq\left(1-\delta^*\right)\|e\|_A^2
$$

we seek a sharp bound

$$
\|S T\|_A^2:=\sup _{e \neq 0} \frac{\|S T e\|_A^2}{\|e\|_A^2}=1-\delta^*
$$

## Sufficient conditions

What should we assume on relaxation and interpolation?

- relaxation is effective on the range of coarse grid correction

There exists $\delta>0$ such that
$$
\|S T \boldsymbol{e}\|_A^2 \leq(1-\delta)\|T \boldsymbol{e}\|_A^2 \text { for all } \boldsymbol{e} .
$$
Then, since $T$ is an $A$-orthogonal projector
$$
\|S T \boldsymbol{e}\|_A^2 \leq(1-\delta)\|\boldsymbol{e}\|_A^2 \text { for all } \boldsymbol{e}
$$

- No side effects on the range of interpolation

$$
\|S \boldsymbol{v}\|_A^2 \leq\|\boldsymbol{v}\|_A^2 \quad \text { for all } \boldsymbol{v} \perp \text{range}(T)
$$

Combine these into an assumption
$$
\|S v\|_A^2 \leq\|S T v\|_A^2+\|S(I-T) v\|_A^2 \leq(1-\delta)\|T v\|_A^2+\|(I-T) v\|_A^2 =\|v\|_A^2-\delta\|T v\|_A^2
$$
**Theorem**
If there exists $\delta>0$ so that
$$
\|S e\|_A^2 \leq\|e\|_A^2-\delta\|T e\|_A^2 \quad \text { for all } e,
$$

then

$$
\|S T\|_A^2 \leq 1-\delta
$$
To be sharp, the largest $\delta$, say $\hat{\delta}$

$$
\hat{\delta}=\inf _{e: T e \neq 0} \frac{\|e\|_A^2-\|S e\|_A^2}{\|T e\|_A^2},
$$

should be $\delta^*$.

Proof

Since $T \boldsymbol{e}=\mathbf{0}$ gives $\|S T \boldsymbol{e}\|_A=0$ :

$$
\|S T\|_A^2=\sup _{e: T e \neq 0} \frac{\|S T e\|_A^2}{\|e\|_A^2}=\sup _{e: T e \neq 0} \frac{\|S T e\|_A^2}{\|T e\|_A^2+\|(I-T) e\|_A^2}
$$

Let $\hat{\boldsymbol{e}}$ be the $\arg\sup$
Then $T \hat{\boldsymbol{e}}$ is also at the supremum.
Thus we have an error at the supremum with $(I-T) \hat{e}=0$
$$
\|S T\|_A^2=\sup _{e: T e \neq 0} \frac{\|S T e\|_A^2}{\|e\|_A^2}=\sup _{e: T e \neq 0} \frac{\|S(T e+(I-T) e)\|_A^2}{\|T e\|_A^2}=\sup _{e: T e \neq 0} \frac{\|S e\|_A^2}{\|T e\|_A^2},
$$


And

$$
1-\|S T\|_A^2=\inf _{e: T e \neq \mathbf{0}} \frac{\|T e\|_A^2-\|S e\|_A^2}{\|T e\|_A^2}=\inf _{e: T e \neq \mathbf{0}} \frac{\|e\|_A^2-\|S e\|_A^2}{\|T e\|_A^2}=\hat{\delta}
$$
The worst $\delta$ is sharp

$$
\hat{\delta}=\inf _{e: T e \neq 0} \frac{\|e\|_A^2-\|S e\|_A^2}{\|T e\|_A^2},
$$
The early theory split this in two part

For some $g(\boldsymbol{e})$ define $\delta, \alpha_g$, and $\beta_g$ as in
$$
\delta(\boldsymbol{e})=\underbrace{\frac{\|\boldsymbol{e}\|_A^2-\|S \boldsymbol{e}\|_A^2}{g(\boldsymbol{e})}}_{\alpha_g(\boldsymbol{e})} \underbrace{\frac{g(\boldsymbol{e})}{\|T \boldsymbol{e}\|_A^2}}_{1 / \beta_g(\boldsymbol{e})}
$$
Consider the smallest $\alpha_g$ and the largest $\beta_g$ :

$$
\hat{\alpha}_g=\inf _{\boldsymbol{e}: g(\boldsymbol{e}) \neq \mathbf{0}} \alpha_g(\boldsymbol{e}) \quad \hat{\beta}_g=\sup _{\boldsymbol{e}: g(\boldsymbol{e}) \neq \mathbf{0}} \beta_g(\boldsymbol{e})
$$
For $\boldsymbol{e}$ such that $g(T \boldsymbol{e}) \neq 0$,

$$
\begin{aligned}
\|S T \boldsymbol{e}\|_A^2 & \leq\|T \boldsymbol{e}\|_A^2-\hat{\alpha}_g g(T \boldsymbol{e}) \leq\|T \boldsymbol{e}\|_A^2-\frac{\hat{\alpha}_g}{\hat{\beta}_g}\|T \boldsymbol{e}\|_A^2=\left(1-\frac{\hat{\alpha}_g}{\hat{\beta}_g}\right)\|T \boldsymbol{e}\|_A^2 \\
& \leq\left(1-\frac{\hat{\alpha}_g}{\hat{\beta}_g}\right)\|\boldsymbol{e}\|_A^2
\end{aligned}
$$
Ok, so this is generally worse than the sharp bound

$$
\|S T\|_A=\sqrt{1-\hat{\delta}} \leq \sqrt{1-\frac{\hat{\alpha}_g}{\hat{\beta}_g}}
$$

( $\alpha_g$ and $\beta_g$ are not generally simultaneously satisfied)

Early works, e.g. Ruge-Stüben 1987, use $g(\boldsymbol{e})=\|\boldsymbol{e}\|_{A D^{-1} A}^2$ Or the weaker form $g(\boldsymbol{e})=\|T \boldsymbol{e}\|_{A D^{-1} A}^2$

## Smoothing and Approximation

if there exists $\bar{\alpha}_g>0$ such that
$$
\|S e\|_A^2 \leq\|e\|_A^2-\bar{\alpha}_g g(e) \quad \text { for all } e \quad \text { (smoothing) }
$$
and there exists $\bar{\beta}_g>0$ such that
$$
\|T \boldsymbol{e}\|_A^2 \leq \bar{\beta}_g g(T \boldsymbol{e}) \quad \text { for all } \boldsymbol{e} \quad \text { (approximation) }
$$
then $\|S T\|_A \leq \sqrt{1-\bar{\alpha}_g / \bar{\beta}_g}$

Select $g(\boldsymbol{e})=\|\boldsymbol{e}\|_{A D^{-1} A}^2$
Since $T$ is an $A$-orthogonal projection we have
$$
\|T e\|_A=\inf _{e_c}\left\|e-P e_c\right\|_A
$$
(strong approximation) Assume there is a $\bar{\beta}_s$ such that
$$
\inf _{\boldsymbol{e}_c}\left\|\boldsymbol{e}-P \boldsymbol{e}_c\right\|_A^2 \leq \bar{\beta}_s\|\boldsymbol{e}\|_{A D^{-1} A}^2 \quad \text { for all } \boldsymbol{e} .
$$
The weaker version looks like (for some $\hat{\beta}$ )

$$
\|T e\|_A^2 \leq \bar{\beta}\|T e\|_{A D^{-1} A} \quad \text { for all } e .
$$


Weaker, means weaker norm. And we can make this a bit more practical. The range of $T$ is $A$-orthogonal to the range of $P$, so

$$
\begin{aligned}
\|T \boldsymbol{e}\|_A^2 & =\langle A T \boldsymbol{e}, T \boldsymbol{e}\rangle=\left\langle A T \boldsymbol{e}, T \boldsymbol{e}-P \boldsymbol{e}_c\right\rangle \\
& \leq\|T \boldsymbol{e}\|_{A D^{-1} A}\left\|T \boldsymbol{e}-P \boldsymbol{e}_c\right\|_D .
\end{aligned}
$$

(weak approximation) Assume that

$$
\inf _{\boldsymbol{e}_c}\left\|\boldsymbol{e}-P \boldsymbol{e}_c\right\|_D^2 \leq \bar{\beta}_w\|\boldsymbol{e}\|_A^2 \quad \text { for all } \boldsymbol{e},
$$
