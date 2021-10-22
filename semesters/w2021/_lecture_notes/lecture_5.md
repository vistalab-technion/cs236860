---
title: "Lecture 5: Statistical estimation"
excerpt: Inverse problems, statistical estimation
author: alex
copyright: alex
---


## Inverse problems and estimation

Let $\mathpzc{F}$ be a latent WSS random signal of which a degraded
observation $\mathpzc{Y} = \mathcal{H} \mathpzc{F} + \mathpzc{N}$ is
given, with $\mathcal{H}$ being a known invertible LSI system (e.g.,
out-of-focus or motion blur), and $\mathpzc{N}$ is zero-mean additive
WSS white noise with $S_ \mathpzc{N}  = \sigma_ \mathpzc{N}^2$
statistically independent of $\mathpzc{F}$ (we will get back to the
question of how to model the noise more realistically). From our
previous derivations, we straightforwardly have
$S_ \mathpzc{Y} = |H|^2 S_ \mathpzc{F} + S_ \mathpzc{N}$ and
$S_ {\mathpzc{F} \mathpzc{Y}} = H^\ast S_ \mathpzc{F}$.

Our goal is to invert the action of $\mathcal{H}$ obtaining the latent
signal $\mathpzc{F}$ from the observation $\mathpzc{Y}$. Such a problem
is known as (non-blind) *deconvolution* and falls into a more general
category of *inverse problems*. In the absence of noise
($\sigma_ \mathpzc{N}=0$), $\mathpzc{F}  = \mathcal{H}^{-1} \mathpzc{Y}$,
where $\mathcal{H}^{-1}$ is the inverse system with the frequency
response $1/H({\bm{\mathrm{\xi}}})$. However, in practice this is very
often a very bad idea. Even if $\mathcal{H}$ is invertible in theory,
its inverse might amplify even the slightest noise to unreasonable
proportions. We will therefore phrase our problem as an *estimation*
problem: find an LSI system $\mathcal{G}$ such that the estimator
$\hat{\mathpzc{F}} = \mathcal{G} \mathpzc{Y}$ is optimally close to the
true $\mathpzc{F}$ in the sense of some error criterion. We define the
*error signal*

$$
\mathpzc{E} = \mathpzc{F} - \hat{\mathpzc{F}} = \mathpzc{F} - \mathcal{G} \mathpzc{Y}
$$

(which is straightforwardly WSS), and seek a system $\mathcal{G}$
satisfying

$$\mathcal{G}_ \ast = \mathrm{arg}\min_ {\mathcal{G}} \epsilon(\mathpzc{E}),$$

where $\epsilon$ is some error criterion. In what follows, we discuss
several choices of $\epsilon$ setting stage to a variety of estimation
frameworks. We also discuss more general inverse problems.

## Wiener filter


A very popular pragmatic choice is the *mean squared error* (MSE),

$$\epsilon = \mathbb{E} \mathpzc{E}({\bm{\mathrm{x}}})^2 = R_ \mathpzc{E}({\bm{\mathrm{0}}}).$$

This choise leads to *minimum MSE* (MMSE) estimators. Using the fact
that 

$$\begin{aligned}
R_ \mathpzc{E}({\bm{\mathrm{x}}}) &= \mathbb{E} \mathpzc{E}({\bm{\mathrm{0}}}) \mathpzc{E}({\bm{\mathrm{x}}}) = \mathbb{E} \left( (\mathpzc{F}({\bm{\mathrm{0}}}) - \hat{\mathpzc{F}}({\bm{\mathrm{0}}}) ) (\mathpzc{F}({\bm{\mathrm{x}}}) - \hat{\mathpzc{F}}({\bm{\mathrm{x}}}) ) \right) \\
&= \mathbb{E} \left( \mathpzc{F}({\bm{\mathrm{0}}})\mathpzc{F}({\bm{\mathrm{x}}}) + \hat{\mathpzc{F}}({\bm{\mathrm{0}}}) \hat{\mathpzc{F}}({\bm{\mathrm{x}}}) - {\mathpzc{F}}({\bm{\mathrm{0}}}) \hat{\mathpzc{F}}({\bm{\mathrm{x}}})  - \hat{\mathpzc{F}}({\bm{\mathrm{0}}}) {\mathpzc{F}}({\bm{\mathrm{x}}})\right)  \\
&= R_ \mathpzc{F}({\bm{\mathrm{x}}}) + R_ {\hat{\mathpzc{F}}}({\bm{\mathrm{x}}}) -  R_ {\mathpzc{F} \hat{\mathpzc{F}}}({\bm{\mathrm{x}}}) -  R_ {\hat{\mathpzc{F}} \mathpzc{F} }({\bm{\mathrm{x}}}), \end{aligned}$$

we can write 

$$\begin{aligned}
\epsilon &=& R_ \mathpzc{E}({\bm{\mathrm{0}}}) = \int_ {\mathbb{R}^d}  S_ \mathpzc{E}({\bm{\mathrm{\xi}}}) \, e^{2\pi \ii\, {\bm{\mathrm{\xi}}}^\Tr {\bm{\mathrm{0}}} } d{\bm{\mathrm{\xi}}} =
 \int_ {\mathbb{R}^d}  ( S_ \mathpzc{F}({\bm{\mathrm{\xi}}}) + S_ {\hat{\mathpzc{F}}}({\bm{\mathrm{\xi}}}) -  S_ {\mathpzc{F} \hat{\mathpzc{F}}}({\bm{\mathrm{\xi}}}) -  S^\ast_ {\mathpzc{F} \hat{\mathpzc{F}}}({\bm{\mathrm{\xi}}}) ) d{\bm{\mathrm{\xi}}}.\end{aligned}$$


Substituting $ S_ \hat{\mathpzc{F}} = |G|^2  S_ {\mathpzc{Y}}$ and
$S_ {\mathpzc{F} \hat{\mathpzc{F}}} = G^\ast S_ {\mathpzc{F} \mathpzc{Y} }$
yields 

$$\begin{aligned}
\epsilon &=& \int_ {\mathbb{R}^d}  ( S_ \mathpzc{F}({\bm{\mathrm{\xi}}}) + |G({\bm{\mathrm{\xi}}})|^2 S_ \mathpzc{Y}({\bm{\mathrm{\xi}}}) - G^\ast({\bm{\mathrm{\xi}}}) S_ {\mathpzc{F}\mathpzc{Y} }({\bm{\mathrm{\xi}}})
- G({\bm{\mathrm{\xi}}}) S^\ast_ {\mathpzc{F}\mathpzc{Y} }({\bm{\mathrm{\xi}}})  ) d\xi.\end{aligned}$$

In order to minimize $\epsilon$ over $G({\bm{\mathrm{\xi}}})$ it is
therefore sufficient to minimize the above integrand for every
${\bm{\mathrm{\xi}}}$ individually; furthermore since the first term in
the integrand does not depend on $G$, we aim at minimizing (for every
${\bm{\mathrm{\xi}}}$, which is omitted for convenience)

$$S_ {\mathpzc{F}} + GG^\ast S_ \mathpzc{Y} - G^\ast S_ {\mathpzc{F}\mathpzc{Y} }  - G S^\ast_ {\mathpzc{F}\mathpzc{Y} } =  GG^\ast S_ \mathpzc{Y} - G^\ast S_ {\mathpzc{F}\mathpzc{Y} }  - G S^\ast_ {\mathpzc{F}\mathpzc{Y} } + \mathrm{const}.$$

(Note that the first term does not depend on $G$).

Observe that at frequencies where $S_ \mathpzc{Y}$ vanishes, it can be
shown that $S_ \mathpzc{FY}$ vanishes as well (Cauchy-Schwarz
inequality). Hence, at those frequencies we may arbitrarily set $G$ to
zero. Otherwise, using the fact that $S_ \mathpzc{Y}$ is real and
non-negative, we can write 

$$\begin{aligned}
GG^\ast S_ \mathpzc{Y} - G^\ast S_ {\mathpzc{F}\mathpzc{Y} }  - G S^\ast_ {\mathpzc{F}\mathpzc{Y} }  =& \left( G S_ \mathpzc{Y}^{\frac{1}{2}} - \frac{ S_ {\mathpzc{F}\mathpzc{Y} } }{S_ \mathpzc{Y}^{\frac{1}{2}}  } \right)\left( G^\ast S_ \mathpzc{Y}^{\frac{1}{2}} - \frac{ S^\ast_ {\mathpzc{F}\mathpzc{Y} } }{S_ \mathpzc{Y}^{\frac{1}{2}}  } \right) - \frac{ S_ \mathpzc{FY}  S_ \mathpzc{FY}^\ast }{ S_ \mathpzc{Y} } \\
=& \left|   G S_ \mathpzc{Y}^{\frac{1}{2}} - \frac{ S_ {\mathpzc{F}\mathpzc{Y} } }{S_ \mathpzc{Y}^{\frac{1}{2}}  }  \right|^2 + \mathrm{const}.\end{aligned}$$


In the absence of other constraints (such as, e.g., bounded spatial
support), the minimizer of the above expression is simply


$$G_ \ast({\bm{\mathrm{\xi}}}) =  \frac{ S_ {\mathpzc{F}\mathpzc{Y} } ({\bm{\mathrm{\xi}}}) }{S_ \mathpzc{Y}({\bm{\mathrm{\xi}}})  }.$$


This result is known as the *Wiener filter* after the mathematician
Norbert Wiener. Note that since both $S_ {\mathpzc{F}\mathpzc{Y} }$ and
$S_ \mathpzc{Y}$ are Fourier transforms or real-valued functions, they
are conjugate symmetric, and so is their ratio. Consequently,
$G_ \ast({\bm{\mathrm{\xi}}})$ is the frequency response of an LSI system
with a real-valued impulse response.

### Estimation MSE

Substituting $G_ \ast$ into the MSE expression yields the error achieved
by the Wiener filter, 

$$\begin{aligned}
\epsilon_ \ast =& \int_ {\mathbb{R}^d}  S_ \mathpzc{E}({\bm{\mathrm{\xi}}})  d{\bm{\mathrm{\xi}}} = \int_ {\mathbb{R}^d}  \left(S_ \mathpzc{F}({\bm{\mathrm{\xi}}})  - G_ \ast({\bm{\mathrm{\xi}}}) S_ \mathpzc{FY}^\ast({\bm{\mathrm{\xi}}}) \right)  d{\bm{\mathrm{\xi}}} \\
=& \int_ {\mathbb{R}^d}  S_ \mathpzc{F}({\bm{\mathrm{\xi}}}) \left( 1 - \frac{ S_ \mathpzc{FY}({\bm{\mathrm{\xi}}}) S_ \mathpzc{FY}^\ast({\bm{\mathrm{\xi}}}) }{ S_ \mathpzc{Y}({\bm{\mathrm{\xi}}}) S_ \mathpzc{F}({\bm{\mathrm{\xi}}}) } \right)  d{\bm{\mathrm{\xi}}} =
\int_ {\mathbb{R}^d}  S_ \mathpzc{F}({\bm{\mathrm{\xi}}}) \left( 1 - \rho({\bm{\mathrm{\xi}}}) \rho^\ast({\bm{\mathrm{\xi}}}) \right)  d{\bm{\mathrm{\xi}}},\end{aligned}$$


where

$$\rho({\bm{\mathrm{\xi}}}) = \frac{ S_ \mathpzc{FY}({\bm{\mathrm{\xi}}})   }{ \sqrt{ S_ \mathpzc{Y}({\bm{\mathrm{\xi}}}) S_ \mathpzc{F}({\bm{\mathrm{\xi}}}) } }$$

plays the role of frequency-wise correlation coefficient.

### Orthogonality

Observe that at every frequency,

$$S_ {\hat{\mathpzc{F}}\mathpzc{Y} } = G_ \ast S_ {\mathpzc{Y} } = \frac{ S_ {\mathpzc{F}\mathpzc{Y} }  }{S_ \mathpzc{Y}   } S_ {\mathpzc{Y} } = S_ {\mathpzc{FY} },$$

from where $R_ {\hat{\mathpzc{F}}\mathpzc{Y} } = R_ {\mathpzc{FY} }$.
Hence,

$$R_ {\mathpzc{E}\mathpzc{Y} }({\bm{\mathrm{\tau}}}) = R_ {(\mathpzc{F}-\hat{\mathpzc{F}})\mathpzc{Y} }({\bm{\mathrm{\tau}}}) = R_ {\mathpzc{F}\mathpzc{Y} }({\bm{\mathrm{\tau}}}) - R_ {\hat{\mathpzc{F}}\mathpzc{Y} }({\bm{\mathrm{\tau}}}) = 0$$


for every ${\bm{\mathrm{\tau}}}$. This can be stated as the following
*orthogonality property*: $\mathpzc{E} \perp \mathpzc{Y}$, that is, the
estimation error is orthonogonal to the data. Orthogonality is the
hallmark of $\ell_ 2$-optimality, which is the case of MMSE estimators.

### Optimal deconvolution

Note that the Wiener filter expression is general and does not assume
any specific relation between the latent signal $\mathpzc{F}$ and the
observation $\mathpzc{Y}$ except the assumption of joint wide-sense
stationarity and that their cross-spectrum is known. In the specific
case of $\mathpzc{Y} = \mathcal{H} \mathpzc{F} + \mathpzc{N}$ with
statistically-independent additive zero-mean noise $\mathpzc{N}$, we
have

$R_ {\mathpzc{Y}}({\bm{\mathrm{\tau}}}) = S_ {(\mathcal{H}\mathpzc{F})}({\bm{\mathrm{\tau}}}) + R_ {\mathpzc{N}}({\bm{\mathrm{\tau}}})$


and can write

$$S_ {\mathpzc{Y}}({\bm{\mathrm{\xi}}}) = |H({\bm{\mathrm{\xi}}})|^2 S_ {\mathpzc{F}}({\bm{\mathrm{\xi}}}) + S_ {\mathpzc{N}}({\bm{\mathrm{\xi}}}).$$


Similarly,

$$S_ {\mathpzc{FY}}({\bm{\mathrm{\xi}}}) = H^\ast({\bm{\mathrm{\xi}}}) S_ {\mathpzc{F}}({\bm{\mathrm{\xi}}}).$$


Substituting into the Wiener filter expression we obtain after
elementary algebraic manipulations

$$G_ \ast({\bm{\mathrm{\xi}}}) =  \frac{  H^\ast({\bm{\mathrm{\xi}}})S_ 
{\mathpzc{F}}({\bm{\mathrm{\xi}}})  }{    H({\bm{\mathrm{\xi}}})  H^\ast({\bm{\mathrm{\xi}}})   S_ {\mathpzc{F}}({\bm{\mathrm{\xi}}}) + S_ {\mathpzc{N}}({\bm{\mathrm{\xi}}})   }.$$

Let us define the *signal-to-noise ratio* (SNR) at frequency
${\bm{\mathrm{\xi}}}$ as

$$\mathrm{SNR}({\bm{\mathrm{\xi}}}) = \frac{ S_ {\mathpzc{F}}({\bm{\mathrm{\xi}}})  }{ S_ {\mathpzc{N}}({\bm{\mathrm{\xi}}})  }.$$


Note that in the particular case of white noise with the spectrum
$S_ \mathpzc{N}  = \sigma_ \mathpzc{N}^2$,

$$\mathrm{SNR}({\bm{\mathrm{\xi}}}) = \frac{ S_ {\mathpzc{F}}({\bm{\mathrm{\xi}}})  }{ \sigma^2_ {\mathpzc{N}}   }.$$

The filter expression can be written as

$$G_ \ast({\bm{\mathrm{\xi}}}) =  \frac{  H^\ast({\bm{\mathrm{\xi}}}) \, \mathrm{SNR}({\bm{\mathrm{\xi}}})  }{    H({\bm{\mathrm{\xi}}})  H^\ast({\bm{\mathrm{\xi}}})  \,  \mathrm{SNR}({\bm{\mathrm{\xi}}}) + 1   }.$$


At frequencies where $\mathrm{SNR}({\bm{\mathrm{\xi}}}) \gg 1$ (signal
is much stronger than the noise),


$$G_ \ast({\bm{\mathrm{\xi}}}) \approx  \frac{  H^\ast({\bm{\mathrm{\xi}}}) \, \mathrm{SNR}({\bm{\mathrm{\xi}}})  }{    H({\bm{\mathrm{\xi}}})  H^\ast({\bm{\mathrm{\xi}}})  \,  \mathrm{SNR}({\bm{\mathrm{\xi}}})  } = \frac{1}{H({\bm{\mathrm{\xi}}})},$$


which is exactly the inverse system. On the other hand, at frequencies
where $\mathrm{SNR}({\bm{\mathrm{\xi}}})$ approaches zero, (noise much
stronger than the signal), the numerator becomes dominant and
$G_ \ast({\bm{\mathrm{\xi}}})$ also approaches zero.

## Maximum likelihood estimators

### Kullback-Leibler divergence

Let $P$ and $Q$ be two probability measures (such that $P$ is absolutely
continuous w.r.t. $Q$). Then, the *Kullback-Leibler divergence* from Q
to P is defined as 

$$D(P || Q) = \int_ {} \, \log \frac{dP}{dQ} \, dP.$$


In other words, it is the expectation of the logarithmic differences
between the probabilities $P$ and $Q$ when the expectation is taken over
$P$. The divergence can be thought of as an (asymmetric) distance
between the two distributions.

### Maximum likelihood

Since $\mathpzc{Y}= \mathcal{H} \mathpzc{F} + \mathpzc{N}$, we can
assert that the distribution of the measurement $\mathpzc{Y}$ given the
latent signal $\mathpzc{F}$ is simply the distribution of $\mathpzc{N}$
at $\mathpzc{N} = \mathpzc{Y}- \mathcal{H} \mathpzc{F}$,

$$P_ {\mathpzc{Y} | \mathpzc{F}}( y | f  ) = P_ {\mathpzc{N}}( y - \mathcal{H} f ).$$


Assuming i.i.d. noise, the latter simplifies to a product of
one-dimensional measures. Note that this is essentially a parametric
family of distributions – each choice of $f$ yields a distribution
$P_ {\mathpzc{Y} | \mathpzc{F} = f}$ of $\mathpzc{Y}$. For the time
being, let us treat the notation $\mathpzc{Y} | \mathpzc{F}$ just as a
funny way of writing.

Given an estimate $\hat{f}$ of the true realization $f$ of
$\mathpzc{F}$, we can measure its “quality” by measuring the distance
from $P_ {\mathpzc{Y} | \mathpzc{F}=\hat{f}}$ to the true distribution
$P_ {\mathpzc{Y} | \mathpzc{F}=f}$ that created $\mathpzc{Y}$, and try to
minimize it. Our estimator of $f$ can therefore be written as


$$\hat{f} = \mathrm{arg}\min_ {\hat{\mathpzc{F}}} D(P_ {\mathpzc{Y} | \mathpzc{F}=f}  ||  P_ {\mathpzc{Y} | \mathpzc{F}=\hat{f}} ),$$


where we used the Kullback-Leibler divergence to quantify the distance
between the distributions. Note that we treat the quantity to be
estimated as a deterministic parameter rather than a stochastic
quantity.

Let us have a closer look at the minimization objective

$$ 
D(P_ {\mathpzc{Y} | \mathpzc{F}=f}  ||  P_ {\mathpzc{Y} | \mathpzc{F}=\hat{f}}  ) = \mathbb{E}_ { \mathpzc{Y} \sim P_ {\mathpzc{Y} | \mathpzc{F}=f}   }   \log\left(  \frac{P_ {\mathpzc{Y} | \mathpzc{F}=f}  }{ P_ {\mathpzc{Y} | \mathpzc{F}=\hat{f} }  } \right) 
=\mathbb{E}_ { \mathpzc{Y} \sim P_ {\mathpzc{Y} | \mathpzc{F} = f}   }   \log P_ {\mathpzc{Y} | \mathpzc{F} = f}  -\mathbb{E}_ { \mathpzc{Y} \sim P_ {\mathpzc{Y} | \mathpzc{F} = f}   }   \log P_ {\mathpzc{Y} | \mathpzc{F}=\hat{f} }.
$$


Note that the first term (that can be recognized as the entropy of
$\log P_ {\mathpzc{Y} | \mathpzc{F}=f}$) does not depend on the
minimization variable; hence, we have


$$\hat{f} = \mathrm{arg}\min_ {\hat{f}} \, \mathbb{E}_ { \mathpzc{Y} \sim P_ {\mathpzc{Y} | \mathpzc{F}=f}   }   \left( - \log P_ {\mathpzc{Y} | \mathpzc{F}=\hat{f}} \right).$$


Let us now assume that $\mathpzc{Y}$ is observed at some set of $N$
spatial locations
$\{ {\bm{\mathrm{x}}}_ 1, \dots, {\bm{\mathrm{x}}}_ N \}$; we will denote
$\mathpzc{Y}_ i = \mathpzc{Y}({\bm{\mathrm{x}}}_ i)$. In this case, we can
use p.d.f.s to write


$$-\frac{1}{N} \log f_ {\mathpzc{Y} | \mathpzc{F}=f} (y_ 1,\dots, y_ N ) = - \frac{1}{N} \sum_ {i=1}^N \log f_ \mathpzc{N} (y_ i - (\mathcal{H} f)({\bm{\mathrm{x}}}_ i) ) = L(y_ 1,\dots,y_ N | f).$$


This function is known as the *negative log likelihood* function. By the
law of large numbers, when $N$ approaches infinity,

$$L(y_ 1,\dots,y_ N | f) \rightarrow \mathbb{E}_ { \mathpzc{Y} \sim P_ {\mathpzc{Y} | \mathpzc{F}=f}   }   
\left( - \log P_ {\mathpzc{Y} | \mathpzc{F}=\hat{f}} \right).$$


Behold our minimization objective!

To recapitulate, recall that we started with minimizing the discrepancy
between the latent parametric distribution that generated the
observation and that associated with our estimator. However, a closer
look at the objective revealed that it is the limit of the negative log
likelihood when the sample size goes to infinity. The minimization of
the Kullback-Leibler divergence is equivalent to maximization of the
likelihood of the data coming from a specific parametric distribution,

$$\hat{f} = \mathrm{arg}\max_ {\hat{f}} \, P(  \mathpzc{Y}=y |  \mathpzc{F}=f ).$$

For this reason, the former estimator is called *maximum likelihood*
(ML).

### ML deconvolution

Back to our deconvolution problem. Assuming white Gaussian distribution
of the noise with zero mean and variance $\sigma^2_ {\mathpzc{N}}$ yields

$$- \log f_ \mathpzc{N}(n) = \mathrm{const} + \frac{n^2 }{2 \sigma^2_ {\mathpzc{N}}};$$


this in turn gives rise to the following negative log likelihood
function


$$L(y_ 1,\dots, y_ N) = \frac{1}{2 N \sigma^2_ {\mathpzc{N} }} \sum_ {i=1}^N (y_ i - (\mathcal{H} f)({\bm{\mathrm{x}}}_ i) )^2.$$


In the limit case, we minimize


$$\int_ {\mathbb{R}^d} (y({\bm{\mathrm{x}}}) - (\mathcal{H} f) ({\bm{\mathrm{x}}}))^2 d{\bm{\mathrm{x}}} = \|  y - \mathcal{H} f \|^2_ {L^2(\mathbb{R}^d)},$$


which by Parseval’s identity is equal to


$$\|  Y - H F \|^2_ {L^2(\mathbb{R}^d)} = \int_ {\mathbb{R}^d} |Y({\bm{\mathrm{\xi}}}) -H({\bm{\mathrm{\xi}}})  F({\bm{\mathrm{\xi}}})|^2 d{\bm{\mathrm{\xi}}}.$$

The latter integral is obviously mimimized by $F = \frac{Y}{H}$, which
we know is a very bad idea in practice.

## Maximum *a posteriori* estimators


### Conditional distributions

Before treating maximum *a posteriori* estimation, we need to briefly
introduce the important notion of conditioning and conditional
distributions. Recall our construction of a probability space comprising
the triplet $\Omega$ (the sample space), $\Sigma$ (the Borel sigma
algebra), and $P$ (the probability measure). Let $X$ be a random
variable and $B \subset \Sigma$ a sub sigma-algebra of $\Sigma$. We can
then define the *conditional expectation of $\mathpzc{X}$ given $B$* as
a random variable $\mathpzc{Z} = \mathbb{E} \mathpzc{X} | B$ satisfying
for every $E \in B$ 

$$\int_ A \mathpzc{Z} dP = \int_ A \mathpzc{X} dP.$$

(we are omitting some technical details such as, e.g., integrability
that $\mathpzc{X}$ has to satisfy).

Given another random variable $\mathpzc{Y}$, we say that it generates a
sigma algebra $\sigma(\mathpzc{Y})$ as the set of pre-images of all
Borel sets in $\mathbb{R}$,

$$\sigma(\mathpzc{Y}) = \{ \mathpzc{Y}^{-1}(A) : A \in \mathbb{B}(\mathbb{R}) \}.$$


We can then use the previous definition to define the conditional
expectation of *$\mathpzc{X}$ given $\mathpzc{Y}$* as


$$\mathbb{E} \mathpzc{X} | \mathpzc{Y} = \mathbb{E} \mathpzc{X} | \sigma(\mathpzc{Y}).$$

Recall that expectation applied to indicator functions can be used to
define probability measures. In fact, for every $E \in \Sigma$, we may
construct the random variable $\ind_ E$, leading to
$P(E) = \mathbb{E} \ind_ E$. We now repeat the same, this time replacing
$\mathbb{E} $ with $\mathbb{E} \cdot | \mathpzc{Y}$. For every
$E \in \Sigma$, $$\varphi(\mathpzc{Y}) = \mathbb{E} \, E | \mathpzc{Y}$$
is a random variable that can be thought of as a transformation of the
random variable $\mathpzc{Y}$ by the function $\varphi$. We denote this
function as $P(E |\mathpzc{Y})$ and refer to it as the (regular)
*conditional probability of $E$ given $\mathpzc{Y}$*. It is easy to show
that for every measurable set $B \subset \mathbb{R}$,


$$\int_ B P(E | \mathpzc{Y}=y) (\mathpzc{Y}_ \ast P)(dy) = P(E \cap \{ \mathpzc{Y} \in B \});$$


Substituting $E = \{ \mathpzc{X} \in B\}$ yields the *conditional
distribution of $X$ given $\mathpzc{Y}$*,


$$P_ {\mathpzc{X} | \mathpzc{Y}} (  B  | \mathpzc{Y}=y) = P(\mathpzc{X} \in B | \mathpzc{Y}=y).$$


It can be easily shown that $P_ {\mathpzc{X} | \mathpzc{Y}}$ is a valid
probability measure on $\Sigma$ and for every pair of measurable sets
$A$ and $B$,


$$\int_ B P_ {\mathpzc{X} | \mathpzc{Y}} (A | \mathpzc{Y}=y) (\mathpzc{Y}_ \ast P)(dy) = P(\{ \mathpzc{X} \in A  \} \cap \{ \mathpzc{Y} \in B \}).$$


If density exists, $P_ {\mathpzc{X} | \mathpzc{Y}}$ can be described
using the *conditional p.d.f.* $f_ {\mathpzc{X} | \mathpzc{Y}}$ and the
latter identity can be rewritten in the form


$$\int_ A \left( \int_ B f_ {\mathpzc{X} | \mathpzc{Y}} (x | y) f_ \mathpzc{Y}(y) dy  \right)  dx = P(\{ \mathpzc{X} \in A  \} \cap \{ \mathpzc{Y} \in B \}) = \int_ A \int_ B f_ {\mathpzc{XY} } (x, y) dxdy.$$


This essentially means that
$f_ {\mathpzc{XY} } (x, y) = f_ {\mathpzc{X} | \mathpzc{Y}} (x | y)  f_ \mathpzc{Y}(y)$.
Integrating w.r.t. $y$ yields the so-called *total probability formula*


$$f_ {\mathpzc{X} } (x) = \int_ \mathbb{R} f_ {\mathpzc{XY} } (x, y) dy = \int_ \mathbb{R} f_ {\mathpzc{X|Y} } (x|y)  f_ \mathpzc{Y}(y) dy.$$


We can also immediately observe that if $\mathpzc{X}$ and $\mathpzc{Y}$
are statistically independent, we have


$$f_ {\mathpzc{XY} } (x, y) = f_ {\mathpzc{X} }(x) f_ {\mathpzc{Y} } (y) =  f_ {\mathpzc{X} | \mathpzc{Y}} (x | y)  f_ \mathpzc{Y}(y),$$


from where $f_ {\mathpzc{X} | \mathpzc{Y}} = f_ {\mathpzc{X}}$. In this
case, conditioning on $\mathpzc{Y}$ does not change our knowledge of
$\mathpzc{X}$.

### Bayes’ theorem

One of the most celebrate (and useful) results related to conditional
distributions is the following theorem named after Thomas Bayes.
Exchanging the roles of $\mathpzc{X}$ and $\mathpzc{Y}$, we have


$$f_ {\mathpzc{XY} }  = f_ {\mathpzc{X} | \mathpzc{Y}}   f_ \mathpzc{Y} = f_ {\mathpzc{Y} | \mathpzc{X}}   f_ \mathpzc{X};$$


re-arranging the terms, we have


$$f_ {\mathpzc{Y} | \mathpzc{X}} = f_ {\mathpzc{X} | \mathpzc{Y}} \, \frac{  f_ \mathpzc{X} }{  f_ \mathpzc{Y} };$$


in terms of probability measures, the equivalent form is


$$P_ {\mathpzc{Y} | \mathpzc{X}} = P_ {\mathpzc{X} | \mathpzc{Y}} \, \frac{  dP_ \mathpzc{X} }{  dP_ \mathpzc{Y} }.$$


### Posterior probability maximization

Recall that in maximum likelihood estimation we treated $\mathpzc{F}$ as
a deterministic parameter and tried to maximize the conditional
probability $P(\mathpzc{Y} | \mathpzc{F})$. Let us now think of
$\mathpzc{F}$ as of a random signal and maximize its probability given
the data,


$$\hat{f}(y) = \mathrm{arg}\max_ { \hat{f} } P_ {\mathpzc{F} | \mathpzc{Y} } ( \mathpzc{F} =  \hat{f} | \mathpzc{Y} = y).$$


Invoking the Bayes theorem yields


$$P_ {\mathpzc{F} | \mathpzc{Y}} = P_ {\mathpzc{Y} | \mathpzc{F} } \, \frac{ dP_ {\mathpzc{F}} }{dP_ {\mathpzc{Y}} }$$


In the Bayesian jargon, $P_ {\mathpzc{F}}$ is called the *prior*
probability, that is, our initial knowledge about $\mathpzc{F}$ before
any observation thereof was obtained; $P_ {\mathpzc{F} | \mathpzc{Y}}$ is
called the *posterior* probability having accounted for the measurement
$\mathpzc{Y}$. Note that the term $P_ {\mathpzc{Y} | \mathpzc{F}}$ is our
good old likelihood. Since we are maximizing the posterior probability,
the former estimator is called *maximum a posteriori* (MAP).

Taking negative logarithm, we obtain


$$-\log P_ {\mathpzc{F} | \mathpzc{Y}} = -\log P_ {\mathpzc{Y} | \mathpzc{F} } -\log P_ \mathpzc{F} +\log P_ {\mathpzc{Y}} = L(\mathpzc{Y} | \mathpzc{F}) - \log P_ {\mathpzc{F}} + \mathrm{const}.$$


This yields the following expression for the MAP estimator


$$\hat{f} = \mathrm{arg}\min_ { \hat{f} } L(\mathpzc{Y} |  \hat{f} ) - \log P_ \mathpzc{F} ( \hat{f} ).$$


The minimization objective looks very similar to what we had in the ML
case; the only difference is that now a *prior* term is added. In the
absence of a good prior, a uniform prior is typically assumed, which
reduces MAP estimation to ML estimation.

### MAP deconvolution

Let us consider for examples our deconvolution problem. As a prior, we
assume that $\mathcal{D} \mathpzc{F}$ is distributed normally i.i.d.
with zero mean and variances $\sigma^2_ \mathrm{D}$, where $\mathcal{D}$
is a derivative of an appropriate order. This leads to the following
objective function


$$\frac{1}{2\sigma^2_ {\mathpzc{N}} }\int_ {\mathbb{R}^d} (y({\bm{\mathrm{x}}}) - (\mathcal{H} f) ({\bm{\mathrm{x}}}))^2 d{\bm{\mathrm{x}}}  + \frac{1}{2\sigma^2_ {D} }\int_ {\mathbb{R}^d} (\mathcal{D} f) ({\bm{\mathrm{x}}})^2 d{\bm{\mathrm{x}}}$$


which by Parseval’s identity is equal to


$$\frac{1}{2\sigma^2_ {\mathpzc{N}} } \int_ {\mathbb{R}^d} |Y({\bm{\mathrm{\xi}}}) -H({\bm{\mathrm{\xi}}})  F({\bm{\mathrm{\xi}}})|^2 d{\bm{\mathrm{\xi}}} + \frac{1}{2\sigma^2_ {D}} \int_ {\mathbb{R}^d} |D({\bm{\mathrm{\xi}}})  F({\bm{\mathrm{\xi}}})|^2 d{\bm{\mathrm{\xi}}}$$


(Note the complex modulus – the integrands in the frequency domain are
complex). The expression is minimized by minimizing at each frequency


$$\begin{aligned}
|Y-HF|^2 + \frac{\sigma^2_ {\mathpzc{N}} }{\sigma^2_ {D}} |DF|^2 &=& Y^\ast Y - Y^\ast HF - Y H^\ast F^\ast + \left( H^\ast H  + \frac{\sigma^2_ {\mathpzc{N}} }{\sigma^2_ {D}} D^\ast D \right) F^\ast F.\end{aligned}$$


Note that this is a convex quadratic function in $F$ having a single
global minimum. Differentiating w.r.t. $F$ and demanding equality to
zero yields, up to the factor of $2$,


$$0 =  -H^\ast Y + \left( H^\ast H  + \frac{\sigma^2_ {\mathpzc{N}} }{\sigma^2_ {D}} D^\ast D \right) F$$


(refer to Chapter 4 in the Matrix Cook Book for the exact treatment of
derivatives w.r.t. a complex variable). The final expression for the
Fourier transform of the estimated signal is given by


$$\hat{F}({\bm{\mathrm{\xi}}}) = \frac{H^\ast({\bm{\mathrm{\xi}}}) }{ H({\bm{\mathrm{\xi}}}) 
H^\ast({\bm{\mathrm{\xi}}}) + \frac{\sigma^2_ {\mathpzc{N}} }{\sigma^2_ {D}} D({\bm{\mathrm{\xi}}}) D^\ast({\bm{\mathrm{\xi}}}) } Y({\bm{\mathrm{\xi}}}).$$


Note the resemblance to the Wiener filter: the ratio
$ \frac{\sigma^2_ {\mathpzc{N}} }{\sigma^2_ {D}}$ controls the amount of
regularization applied to the inverse of $H$. When the noise is strong,
the term $D^\ast D$ dominates; for vanishing noise variance, the
expression reduces to $H^{-1}$. However, the MAP formalizm immediately
allows incorporating more meaningful and realistic priors instead of our
toy example. While not leading to closed-form solutions, they will still
form well-defined minimization problems that can be solved (or at least
attempted to be solved) iteratively. An important philosophical question
arises at this point: What is better – to use an incorrect model that we
know how to solve (e.g., the Wiener filter), or try to develop a better
model that we can hope to solve (e.g., some non-convex MAP estimation
problem)? Practice of the past two decades shows that the second
approach may literally lead (and has, indeed, led) to revolutions.
Without too much exaggeration, we can state that over $99\%$ of recent
studies in image processing essentially revolve around finding a better
expression for the prior probability $P_ {\mathpzc{F} }$.

### Bayesian estimators

Note that the posterior distribution $P_ {\mathpzc{F} | \mathpzc{Y}}$ can
be used to define conditional expectations. Let $\ell( f, \hat{f})$ be
some *loss function* determining how far our estimate $\hat{f}$ is from
the true signal $f$. Then, we can define an error criterion

$$\epsilon = \mathbb{E} \, \ell( \mathpzc{F}, \hat{\mathpzc{F}}) | \mathpzc{Y}$$

and minimize it over estimators of the form $\hat{f} = \hat{f}(y)$.
(Note that the expectation is evaluated w.r.t. the posterior
distribution $P_ {\mathpzc{F} | \mathpzc{Y}}$). This yields to the
so-called *Bayesian estimators*

$$
\hat{f} = \mathrm{arg}\min_ {\hat{f}} \mathbb{E} \, \ell( \mathpzc{F}, \hat{f}(\mathpzc{Y}) ) | \mathpzc{Y}.
$$

The Wiener filter (MMSE estimator) is a particular case of Bayesian
estimators where the loss is set to be Euclidean,
$\ell(f, \hat{f}) = \| f - \hat{f} \|^2_ {L^2(\mathbb{R}^d}$, and
estimators have the form $\hat{f}(y) = h\ast y$.

For the particular choice of the loss function

$$
\ell(f, \hat{f}  = \left\{ \begin{array}{cl}  0 & : | f- \hat{ f } | \le c \\  1 & : \mathrm{else} \end{array} \right.
$$

a $c$ approaches $0$, the Bayesian estimator estimates the mode of the
posterior distribution and thus approaches the MAP estimator, provided
that the distribution of $\mathpzc{F}$ is unimodal.
