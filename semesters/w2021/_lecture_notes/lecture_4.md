---
title: "Lecture 4: Random signals"
excerpt: Elementary probabilty, Random signals
author: alex
copyright: alex
---

## Elementary probability


### Probability measure

We start with a few elementary (and simplified) definitions from the
theory of probability. Let us fix a *sample space* $\Omega = [0,1]$. A
*Borel set* on $\Omega$ is a set that can be formed from open intervals
of the form $(a,b), 0 \le a<b \le 1$, through the operations of
countable union, countable intersection, and set difference. We will
denote the collection of all Borel sets in $\Omega$ as $\Sigma$. It is
pretty straightforward to show that $\Sigma$ contains the empty set, is
closed under complement, and is closed under countable union. Such a set
is known as *$\sigma$-algebra* and its elements (subsets of
$\mathbb{R}$) are referred to as *events*.

A *probability measure* $P$ on $\Sigma$ is a function
$P : \Sigma \rightarrow [0,1]$ satisfying $P(\emptyset) = 0$,
$P(\mathbb{R}) = 1$ and additivity for every countable collection
$\{ E _n \in \Sigma \}$,

$$P\left(  \bigcup _n E _n \right) = \sum _{n} P(E _n).$$

### Random variables

A *random variable* $\mathpzc{X}$ is a *measurable* map
$\mathpzc{X} : \Omega \rightarrow \mathbb{R}$, i.e., a function such
that for every $a$,
$\{ \mathpzc{X} \le a \} = \{ \alpha : \mathpzc{X}(\alpha) \le a  \} \in \Sigma$.
The map $\mathpzc{X}$ pushes forward the probability measure $P$; the
*pushforward measure* $\mathpzc{X} _\ast P$ is given by

$$(\mathpzc{X} _\ast P)(A) = P(\mathpzc{X}^{-1}(A)),$$

where
$\mathpzc{X}^{-1}(A) = \{ \alpha  : X(\alpha) \in A \}$ is the preimage
of $A \subseteq \mathbb{R}$. (In short, we can write
$\mathpzc{X} _\ast P = P\mathpzc{X}^{-1}$). This pushforward probability
measure $\mathpzc{X} _\ast P$ is usually referred to as the *probability
distribution* (or the *law*) of $\mathpzc{X}$.

When the range of $\mathpzc{X}$ is finite or countably infinite, the
random variable is called *discrete* and its distribution can be
described by the *probability mass function* (PMF):

$$f _{\mathpzc{X}}(x) = P(\mathpzc{X}=x),$$ 

which is a shorthand for $P(  \{\alpha : \mathpzc{X}(\alpha) = x \} )$. 
Otherwise, $\mathpzc{X}$ is called a *continuous* random variable.
Any random variable can be described by the *cumulative distribution function* (CDF)

$$F _{\mathpzc{X}}(x) = P({\mathpzc{X} \le x}),$$ 

which is a shorthand for
$F _{\mathpzc{X}}(x) = P(  \{\alpha : \mathpzc{X}(\alpha) \le x \} )$. If
$X$ is absolutely continuous, the CDF can be described by the integral

$$F _{\mathpzc{X}}(x) = \int _{-\infty}^x f _{\mathpzc{X}}(x') dx',$$ 

where the integrand $f _{\mathpzc{X}}$ is known as the *probability density
function* (PDF)[^1].

### Uniform distribution and uniformization

A random variable $\mathpzc(U)$ is said to be *uniformly distributed* on
$[0,1]$ (denoted as $\mathpzc{U} \sim \mathcal{U}[0,1]$) if

$$P(\mathpzc{U} \in [a,b]) = b-a = \lambda([a,b]).$$ 

In other words, the map $\mathpzc{U}$ pushes forward the standard Lebesgue measure on
$[0,1]$, $\mathpzc{U} _\ast P = \lambda$. The corresponding CDF is
$F _\mathpzc{U}(u) = \max\{ 0, \min\{ 1, u \} \}$. Let $\mathpzc{X}$ be
some other random variable characterized by the CDF $F _\mathpzc{X}$. We
define $\mathpzc{U} = F _\mathpzc{X}(\mathpzc{X})$. Let us pick an
arbitrary $x \in \mathbb{R}$ and let $u = F _\mathpzc{X}(x) \in [0,1]$.
From monotonicity of the CDF, it follows that $\mathpzc{U} \le u$ if and
only if $\mathpzc{X} \le x$. Hence,
$F _\mathpzc{U}(u) = P(\mathpzc{U} \le u) = P(\mathpzc{X} \le x) = F _\mathpzc{X}(x) = u$.
We conclude that by transforming a random variable with its own CDF
uniformizes it on the interval $[0,1]$.

Applying the relation in inverse direction, let
$\mathpzc{U} \sim \mathcal{U}[0,1]$ and let $F$ be a valid CDF. Then,
the random variable $\mathpzc{X} = F^{-1}(\mathpzc{U})$ is distributed
with the CDF $F _\mathpzc{U} = F$.

### Expectation

The *expected value* (a.k.a. the *expectation* or *mean*) of a random
variable $\mathpzc{X}$ is given by

$$\mathbb{E} \mathpzc{X} = \int _{\mathbb{R}} \mathrm{id}\, d(\mathpzc{X} _\ast P) =  \int _{\Omega} \mathpzc{X}(\alpha) d\alpha,$$

where the integral is the Lebesgue integral w.r.t. the measure $P$;
whenever a probability density function exists, the latter can be
written as

$$\mathbb{E} \mathpzc{X} = \int _{\mathbb{R}} x  f _{\mathpzc{X}}(x) dx.$$

Note that due to the linearity of integration, the expectation operator
$\mathbb{E}$ is linear. Using the Lebesgue integral notation, we can
write for $E \in \Sigma$

$$P(E) = \int _E dP = \int _\mathbb{R} \ind _E \, dP = \mathbb{E}  \ind _E,$$

where 

$$\ind _E(\alpha) = \left\{
 \begin{array}{ccc}
 1 & : & \alpha \in E \\
 0 & : & \mathrm{otherwise}
 \end{array}
   \right.$$ 

is the *indicator function* of $E$, which is by itself a
random variable. This relates the expectation of the indicator of an
event to its probability.

### Moments

For any measurable function $g : \mathbb{R} \rightarrow \mathbb{R}$,
$\mathpzc{Z} = g(\mathpzc{X})$ is also a random variable with the
expectation

$$\mathbb{E} \mathpzc{Z} = \mathbb{E} g(\mathpzc{X}) =  \int _{\mathbb{R}} g\, dP =  \int _{\mathbb{R}} g(x) f _{\mathpzc{X}}(x) dx.$$

Such an expectation is called a *moment* of $\mathpzc{X}$. Particularly,
the $k$-th order moment is obtained by setting $g(x) = x^k$,

$$\mu _{k}(\mathpzc{X}) = \mathbb{E} \mathpzc{X}^k.$$ 

The expected value itself is the first-order moment of $\mathpzc{X}$, which is often
denoted simply as $\mu _\mathpzc{X} = \mu _{1}(\mathpzc{X})$. The
*central* $k$-th order moment is obtained by setting
$g(x) = (x - \mu _\mathpzc{X})^k$,

$$m _{k}(\mathpzc{X}) = \mathbb{E} ( \mathpzc{X}  - \mathbb{E}  \mathpzc{X})^k.$$

A particularly important central second-order moment is the *variance*

$$\sigma _\mathpzc{X}^2 = \mathrm{Var}\, \mathpzc{X} = m _2 = \mathbb{E} ( \mathpzc{X}  - \mathbb{E}  \mathpzc{X})^2 = \mu _2  ( \mathpzc{X} ) - \mu^2 _\mathpzc{X}.$$

### Joint distribution

A vector $\mathpzcb{X} = (\mathpzc{X} _1, \dots, \mathpzc{X} _n)$ of
random variables is called a *random vector*. Its probability
distribution is defined as before as the pushforward measure
$P = \mathpzcb{X} _\ast \lambda$ Its is customary to treat $\mathpzcb{X}$
as a collection of $n$ random variables and define their *joint CDF* as

$$F _{\mathpzcb{X}}(\bb{x}) = P({\mathpzcb{X} \le \bb{x}}) = P(\mathpzc{X} _1 \le x _1, \dots, \mathpzc{X} _n \le x _n) 
= P(\{ \mathpzc{X} _1 \le x _1 \} \times  \dots \times \{ \mathpzc{X} _n \le x _n \}).$$

As before, whenever the following holds

$$F _{\mathpzcb{X}}(\bb{x}) = \int _{-\infty}^{x _1} \cdots \int _{-\infty}^{x _n}  f _{\mathpzcb{X}}(x _1',\dots, x _n') dx' _1 \cdots dx' _n,$$

the integrand $f _{\mathpzcb{X}}$ is called the *joint PDF* of
$\mathpzcb{X}$. The more rigorous definition as the Radon-Nikodym
derivative

$$f _{\mathpzcb{X}} = \frac{d(\mathpzc{X} _\ast P) }{ d\lambda}$$ 

stays unaltered, only that now $\lambda$ is the $n$-dimensional Lebesgue
measure.

### Marginal distribution

Note that the joint CDF of the sub-vector
$(\mathpzc{X} _2, \dots, \mathpzc{X} _n)$ is given by 

$$\begin{aligned}
F _{\mathpzc{X} _2, \cdots, \mathpzc{X} _n } (x _2,\dots,x _n) =& P(\mathpzc{X} _2 \le x _2, \dots, \mathpzc{X} _n \le x _n) 
= P(\mathpzc{X} _1 \le \infty, \mathpzc{X} _2 \le x _2, \dots, \mathpzc{X} _n \le x _n)  \\
=& F _{\mathpzc{X} _1, \cdots, \mathpzc{X} _n } (\infty, x _2,\dots,x _n).\end{aligned}$$

Such a distribiution is called *marginal* w.r.t. $\mathpzc{X} _1$ and the
process of obtaining it by substituting $x _1 = \infty$ into
$F _{\mathpzcb{X}}$ is called *marginalization*. The corresponding action
in terms of the PDF consists of integration over $x _1$,

$$\begin{aligned}
f _{\mathpzc{X} _2, \cdots, \mathpzc{X} _n } (x _2,\dots,x _n) =& 
\int _\mathbb{R} f _{\mathpzc{X} _1, \cdots, \mathpzc{X} _n } (x _1, x _2,\dots,x _n) dx _1.\end{aligned}$$

### Statistical independence

A set $\mathpzc{X} _1, \dots, \mathpzc{X} _n$ of random variables is
called *statistically independent* if their joint CDF is
coordinate-separable, i.e., can be written as the following tensor
product

$$F _{\mathpzc{X} _1, \cdots, \mathpzc{X} _n} = F _{\mathpzc{X} _1} \otimes \cdots \otimes F _{\mathpzc{X} _n}.$$

An alternative definion can be given in terms of the PDF (whenever it
exists):

$$f _{\mathpzc{X} _1, \cdots, \mathpzc{X} _n} = f _{\mathpzc{X} _1} \otimes \cdots \otimes f _{\mathpzc{X} _n}.$$

We will see a few additional alternative definitions in the sequel.

### Convolution theorem

Let $\mathpzc{X}$ and $\mathpzc{Y}$ be statistically-independent random
variables with a PDF and let $\mathpzc{Z} = \mathpzc{X}+\mathpzc{Y}$.
Then, 

$$\begin{aligned}
F _\mathpzc{Z}(z) =& P(\mathpzc{Z} \le z) = P(X+Y \le z) = \int _{\mathbb{R}} \int _{\infty}^{z-y} f _{\mathpzc{X}\mathpzc{Y}}(x,y) dxdy \\
=& \int _{\mathbb{R}} \int _{\infty}^{z} f _{\mathpzc{X}\mathpzc{Y}}(x'-y,y) dx' dy,\end{aligned}$$


where we changed the variable $x$ to $x' = x+y$. Differentiating w.r.t.
$z$ yields

$$\begin{aligned}
f _\mathpzc{Z}(z)  =& \frac{dF _\mathpzc{Z}(z)}{dz} = \int _{\mathbb{R}} \frac{\partial}{\partial z} \int _{\infty}^{z} f _{\mathpzc{X}\mathpzc{Y}}(x'-y,y) dx' dy = \int _{\mathbb{R}}  f _{\mathpzc{X}\mathpzc{Y}}(z-y,y)  dy.\end{aligned}$$

Since $\mathpzc{X}$ and $\mathpzc{Y}$ are statistically-independent, we
can substitute
$f _{\mathpzc{X}\mathpzc{Y}} = f _{\mathpzc{X}} \otimes f _{\mathpzc{Y}}$
yielding 

$$\begin{aligned}
f _\mathpzc{Z}(z)  =& \int _{\mathbb{R}}  f _{\mathpzc{X}}(z-y) f _{\mathpzc{Y}}(y)  dy = (f _{\mathpzc{X}} \ast f _{\mathpzc{Y}} )(z).\end{aligned}$$

### Limit theorems

Given independent identically distributed (i.i.d.) variables
$\mathpzc{X} _1, \dots, \mathpzc{X} _n$ with mean $\mu$ and variance
$\sigma^2$, we define their *sample average* as

$$\mathpzc{S} _n = \frac{1}{n}( \mathpzc{X} _1 + \cdots + \mathpzc{X} _n ).$$

Note that $\mathpzc{S} _n$ is also a random variable with
$\mu _{\mathpzc{S} _n} = \mu$ and
$\displaystyle{\sigma^2 _{\mathpzc{S} _n} = \frac{\sigma^2}{n}}$. It is
straightforward to see that the variance decays to zero in the limit
$n \rightarrow \infty$, meaning that $\mathpzc{S} _n$ approaches a
deterministic variable $\mathpzc{S} = \mu$. However, a much stronger
result exists: the (strong) *law of large numbers* states that in the
limit $n \rightarrow \infty$, the sample average converges in
probability to the expected value, i.e.,

$$P\left(  \lim _{n \rightarrow \infty} \mathpzc{S} _n = \mu  \right) = 1.$$

This fact is often denoted as
$\mathpzc{S} _n \mathop{\rightarrow}^P \mu$. Furthermore, defining the
normalized deviation from the limit
$\mathpzc{D} _n = \sqrt{n}(\mathpzc{S} _n - \mu)$, the *central limit
theorem* states that $\mathpzc{D} _n$ converges in distribution to
$\mathcal{N}(0,\sigma^2)$, that is, its CDF converges pointwise to that
of the normal distribution. This is often denoted as
$\mathpzc{D} _n \mathop{\rightarrow}^D \mathcal{N}(0,\sigma^2)$.

A slightly more general result is known as the *delta method* in
statistics: if $g :  \mathbb{R} \rightarrow \mathbb{R}$ is a
$\mathcal{C}^1$ function with non-vanishing derivative, then by the
Taylor theorem,

$$g(\mathpzc{S} _n) = g(\mu) + g'(\nu)(\mathpzc{S} _n-\mu) + \mathcal{O}(| \mathpzc{S} _n-\mu |^2),$$

where $\nu$ lies between $\mathpzc{S} _n$ and $\mu$. Since by the law of
large numbers $\mathpzc{S} _n \mathop{\rightarrow}^P \mu$, we also have
$\nu \mathop{\rightarrow}^P \mu$; since $g'$ is continuous,
$g'(\nu) \mathop{\rightarrow}^P g'(\mu)$. Rearranging the terms and
multiplying by $\sqrt{n}$ yields

$$\sqrt{n}( g(\mathpzc{S} _n) - g(\mu) ) = g'(\nu) \sqrt{n}( \mathpzc{S} _n) - \mu ) = g'(\nu) \mathpzc{D} _n,$$

from where (formally, by invoking the Slutsky theorem):

$$\sqrt{n}( g(\mathpzc{S} _n) - g(\mu) ) \mathop{\rightarrow}^D \mathcal{N}(0,g^{\prime} (\mu)^2 \sigma^2).$$

### Joint moments

Given a measurable function
$\bb{g} : \mathbb{R}^n \rightarrow \mathbb{R}^m$, a (joint) moment of a
random vector $\mathpzcb{X} = (\mathpzc{X} _1, \dots, \mathpzc{X} _n)$ is

$$\mathbb{E} \bb{g}(\mathpzcb{X}) = \int \bb{g}(\bb{x}) dP = 
\left(\begin{array}{c} 
 \int g _1(\bb{x}) dP \\
\vdots \\
 \int g _m(\bb{x}) dP
\end{array}
\right)
=
\left(\begin{array}{c} 
\int _{\mathbb{R}^n} g _1(\bb{x}) f _{\mathpzcb{X}}(\bb{x}) d\bb{x}  \\
\vdots \\
\int _{\mathbb{R}^n} g _m(\bb{x}) f _{\mathpzcb{X}}(\bb{x}) d\bb{x}
\end{array}
\right);$$ 

the last term migh be undefined if the PDF does not exist.
The mean of a random vector is simply
$\bb{\mu} _\mathpzcb{X}  = \mathbb{E} \mathpzcb{X}$. Of particular
importance are the second-order joint moments of pairs of random
variables,

$$r _{\mathpzc{X}\mathpzc{Y}} = \mathbb{E} \mathpzc{X}\mathpzc{Y}$$ 

and its central version

$$\sigma^2 _{\mathpzc{X}\mathpzc{Y}} = \mathrm{Cov}(\mathpzc{X},\mathpzc{Y}) = \mathbb{E} \left( (\mathpzc{X} - \mathbb{E} \mathpzc{X} )(\mathpzc{Y}  - \mathbb{E} \mathpzc{Y}) \right) = r _{\mathpzc{X}\mathpzc{Y}} - \mu _\mathpzc{X} \mu _\mathpzc{Y}.$$

The latter quantity is known as the *covariance* of $\mathpzc{X}$ and
$\mathpzc{Y}$.

Two random variables $\mathpzc{X}$ and $\mathpzc{Y}$ with
$r _{\mathpzc{X}\mathpzc{Y}} = 0$ are called *orthogonal*[^2]

The variables with $\sigma^2 _{\mathpzc{X}\mathpzc{Y}} = 0$ are called
*uncorrelated*. Note that for a statistically independent pair
$(\mathpzc{X},\mathpzc{Y})$, 

$$\begin{aligned}
\sigma^2 _{\mathpzc{X}\mathpzc{Y}} =& \int _{\mathbb{R}^2} (x-\mu _\mathpzc{X}) (y-\mu _\mathpzc{Y}) d((\mathpzc{X} \times \mathpzc{Y}) _\ast P) = \int _{\mathbb{R}} (x-\mu _\mathpzc{X}) d(\mathpzc{X} _\ast P) \, \int _{\mathbb{R}} (y-\mu _\mathpzc{Y}) d(\mathpzc{Y} _\ast P) \\
=& \mathbb{E} (\mathpzc{X} - \mathbb{E} \mathpzc{X} ) \cdot \mathbb{E} (\mathpzc{Y}  - \mathbb{E} \mathpzc{Y}) = 0.\end{aligned}$$

However, the converse is not true, i.e., lack of correlation does not
generally imply statistical independence (with the notable exception of
normal variables). If $\mathpzc{X}$ and $\mathpzc{Y}$ are uncorrelated
and furthermore one of them is zero-mean, then they are also orthogonal
(and the other way around).

In general, the *correlation matrix* of a random vector
$\mathpzcb{X} = (\mathpzc{X} _1, \dots, \mathpzc{X} _n)$ is given by

$$\bb{R} _{\mathpzcb{X}} = \mathbb{E}  \mathpzcb{X} \mathpzcb{X}^\Tr;$$

its $(i,j)$-th element is
$(\bb{R} _{\mathpzcb{X}}) _{ij} = \mathbb{E} \mathpzc{X} _i \mathpzc{X} _j$.
Similarly, the *covariance matrix* is defined as the central counterpart
of the above moment,

$$\bb{C} _{\mathpzcb{X}} = \mathbb{E}  (\mathpzcb{X} - \bb{\mu} _\mathpzcb{X} ) (\mathpzcb{X} - \bb{\mu} _\mathpzcb{X} )^\Tr;$$

its $(i,j)$-th element is
$(\bb{C} _{\mathpzcb{X}}) _{ij} =\mathrm{Cov}( \mathpzc{X} _i , \mathpzc{X} _j)$.
Given another random vector
$\mathpzcb{Y} = (\mathpzc{Y} _1, \dots, \mathpzc{Y} _m)$, the
*cross-correlation* and *cross-covariance* matrices are defined as
$\bb{R} _{\mathpzcb{X}\mathpzcb{Y}} = \mathbb{E}  \mathpzcb{X} \mathpzcb{Y}^\Tr$
and
$\bb{C} _{\mathpzcb{X}\mathpzcb{Y}} = \mathbb{E}  (\mathpzcb{X} - \bb{\mu} _\mathpzcb{X} ) (\mathpzcb{Y} - \bb{\mu} _\mathpzcb{Y} )^\Tr$,
respectively.

### Linear transformations

Let $\mathpzcb{X} = (\mathpzc{X} _1, \dots, \mathpzc{X} _n)$ be an
$n$-dimensional random vector, $\bb{A}$ and $m \times n$ deterministic
matrix, and $\bb{b}$ and $m$-dimensional deterministic vector. We define
a random vector $\mathpzcb{Y} = \bb{A} \mathpzcb{X} + \bb{b} $ as the
affine transformation of $\mathpzcb{X}$. Using linearity of the
expectation operator, it is straightforward to show that

$$\begin{aligned}
\bb{\mu} _\mathpzcb{Y}  =& \mathbb{E}(\bb{A} \mathpzcb{X} + \bb{b}) = \bb{A} \bb{\mu} _\mathpzcb{X} + \bb{b} \\
\bb{C} _\mathpzcb{Y}  =& \mathbb{E}(\bb{A} \mathpzcb{X} - \bb{A} \bb{\mu} _\mathpzcb{X} ) (\bb{A} \mathpzcb{X} - \bb{A} \bb{\mu} _\mathpzcb{X} )^\Tr = \bb{A} \bb{C} _\mathpzcb{X} \bb{A}^\Tr \\
\bb{C} _{\mathpzcb{X} \mathpzcb{Y}}  =& \mathbb{E}(\mathpzcb{X} - \bb{\mu} _\mathpzcb{X} ) (\bb{A} \mathpzcb{X} - \bb{A} \bb{\mu} _\mathpzcb{X} )^\Tr  = \bb{C} _\mathpzcb{X} \bb{A}^\Tr.\end{aligned}$$

## Random signals


A *random signal* or a *stochastic process* is a collection of random
variables $\{ \mathpzc{F}(\bb{x}) : \bb{x} \in D \}$ indexed by some set
$D$ (called the domain) and assuming values in some space $S$ (called
the state space). For our purpose, we will assume $D$ to be either
$\mathbb{R}^d$ (in this case we will refer to the process as to a
continuous-domain stochastic process) or $\mathbb{Z}^d$ (discrete-domain
process); the state space $D$ will be assumed either $\mathbb{R}$
(continuous state) or $\mathbb{Z}$ (discrete state). Informally, setting
$D=\mathbb{R}^2$ and $S=\mathbb{R}$ we can think of $\mathpzc{F}$ as of
a random image. When $d>1$, it is customary to call the stochastic
process a *random field*. We will henceforth use the term “random
signal”. In what follows, we will almost always tacitly assume $D$ to be
$\mathbb{R}^d$; the very same reasoning applies to discrete-domain
signals mutatis mutandis.

Formally, a random signal $\mathpzc{F}$ is a function
$\mathpzc{F} : D \times \Omega \rightarrow S$, where $\Omega$ denotes
the sample space. The first argument $\bb{x} \in D$ in
$\mathpzc{F}(\bb{x}, \alpha)$ sets the location in the domain, while the
second argument $\alpha \in \Omega$ is responsible for the randomness.
Fixing some $\alpha \in \Omega$, the resulting deterministic function
$f(\bb{x}) = \mathpzc{F}(\bb{x}, \alpha)$ is called a *realization* or a
*sample function* of $\mathpzc{F}$. Fixing some $\bb{x} \in D$, we
obtain a random variable
$\mathpzc{F}(\alpha) =  \mathpzc{F}(\bb{x}, \alpha)$ describing the
randomness of the signal sampled at a fixed location $\bb{x}$. Note that
for a singletone $D = \{1\}$, the signal is just a random variable,
while for a finite $D = \{1,2,\dots,n\}$ it is a random vector. Random
signals can be therefore thought of as a generalization of random
vectors to “random functions”. However, in such an
infinitely-dimensional case it is not easy to define the standard
notions such as density. Instead, we are going to sample the signal at a
set of points $\{ \bb{x} _1,\dots,\bb{x} _n \} \subset D$ and describe the
distribution of the random vector
$\mathpzcb{F} = (\mathpzc{F}(\bb{x} _1),\dots, \mathpzc{F}(\bb{x} _n))$.
We will define the *finite-dimensional CDF* as

$$F _{\bb{x} _1,\dots,\bb{x} _n}(f _1,\dots,f _n) = P(\mathpzc{F}(\bb{x} _1) \le f _1, \dots, \mathpzc{F}(\bb{x} _n) \le f _n).$$

### Stationarity

A random signal is said (strict sense) *stationary* (SSS), if its
probability distribution is shift-invariant. In other words, if for any
$n$, any $\{ \bb{x} _1,\dots,\bb{x} _n \} \subset D$, and any
$\bb{p} \in D$,
$(\mathpzc{F}(\bb{x} _1+\bb{p}),\dots, \mathpzc{F}(\bb{x} _n+\bb{p}))$ has
the same distribution as
$(\mathpzc{F}(\bb{x} _1),\dots, \mathpzc{F}(\bb{x} _n))$, then
$\mathpzc{F}$ is SSS.

### Auto-correlation and cross-correlation

Moments of random signals can be defined by considering random vectors
obtained from sampling the signal at some set of locations. Of
particular importance will be the following first- and second-order
moments: the *mean* function

$$\mu _{\mathpzc{F}}(\bb{x}) = \mathbb{E} \mathpzc{F}(\bb{x}))$$ 

and the *auto-correlation* function

$$R _{\mathpzc{F}}(\bb{x} _1,\bb{x} _2) = \mathbb{E} \mathpzc{F}(\bb{x} _1) \mathpzc{F}(\bb{x} _2).$$

Given two random signals $\mathpzc{F}$ and $\mathpzc{G}$, we can
similarly define their *cross-correlation* as

$$R _{\mathpzc{F}\mathpzc{G}}(\bb{x} _1,\bb{x} _2) = \mathbb{E} \mathpzc{F}(\bb{x} _1) \mathpzc{G}(\bb{x} _2).$$

It follows from definition that
$R _{\mathpzc{F}}(\bb{x} _1,\bb{x} _2)  = R _{\mathpzc{F}}(\bb{x} _2,\bb{x} _1)$
and
$R _{\mathpzc{G}\mathpzc{F}}(\bb{x} _1,\bb{x} _2)  = R _{\mathpzc{F}\mathpzc{G}}(\bb{x} _2,\bb{x} _1)$.
Two random signals are said to be *orthogonal* if their covariance
function vanishes at every point.

### Wide-sense stationarity

A random signal $\mathpzc{F}$ is called *wide-sense stationary* (WSS) if
its mean and auto-correlation functions are shift-invariant, i.e., for
every $\bb{p} \in D$
$\mu _{\mathpzc{F}}(\bb{x} + \bb{p}) = \mu _{\mathpzc{F}}(\bb{x})$ for
every $\bb{x} \in D$, and
$R _{\mathpzc{F}}(\bb{x} _1 + \bb{p},\bb{x} _2 + \bb{p}) = R _{\mathpzc{F}}(\bb{x} _1,\bb{x} _2)$
for every $\bb{x} _1,\bb{x} _2 \in D$. These conditions immediately
translate to demanding $\mu _{\mathpzc{F}} = \mathrm{const}$ and
$R _{\mathpzc{F}}(\bb{x} _1,\bb{x} _2) = R _{\mathpzc{F}}(\bb{x} _1 - \bb{x} _2)$.
Two WSS random signals $\mathpzc{F}$ and $\mathpzc{G}$ are called
*jointly* WSS if their cross-correlation is shift-invariant, i.e.,
$R _{\mathpzc{F}\mathpzc{G}}(\bb{x} _1,\bb{x} _2) = R _{\mathpzc{F}\mathpzc{G}}(\bb{x} _1 - \bb{x} _2)$.

### Power spectrum density

We start our discussion with the more familiar domain of deterministic
signals. In the signal processing jargon, the *energy* of a
deterministic signal $f$ is defined as

$$E = \| f \|^2 _{L^2(\mathbb{R}^d)} = \int _{\mathbb{R}^d} | f(\bb{x})|^2 d\bb{x}.$$

Due to Parseval’s identity,

$$E = \| \mathcal{F} f \|^2 _{L^2(\mathbb{R}^d)} = \int _{\mathbb{R}^d} | F(\bb{\xi})|^2 d\bb{\xi}$$

(the same is true mutatis mutandis for discrete-domain signals). We can
therefore think of $| F(\bb{\xi})|^2$ as of the *energy density* of $f$
per unit of frequency.

When the signal has infinite energy, we can still define its *average
power* by windowing the signal and normalizing its energy within the
window by the volume of the latter,

$$W = \lim _{T \rightarrow \infty} \frac{1}{T^d} \int _{[-\frac{T}{2},\frac{T}{2}]^d} | f(\bb{x})|^2 d\bb{x} =\lim _{T \rightarrow \infty} \left\|  \frac{1}{T^\frac{d}{2}} \, \mathrm{rect} _T \cdot f \right\|^2 _{L^2(\mathbb{R}^d)}$$

where $\mathrm{rect} _T(\bb{x}) = \mathrm{rect}(\bb{x}/T)$. Defining the
windowed Fourier transform

$$F _T(\bb{\xi}) =  \frac{1}{T^\frac{d}{2}} \, \mathcal{F}\left( \mathrm{rect} _T \cdot f \right)$$

and invoking Parseval’s identity, we have

$$W = \lim _{T \rightarrow \infty} \| F _T(\bb{\xi}) \|^2 _{L^2(\mathbb{R}^d)} = 
\lim _{T \rightarrow \infty}   \int _{\mathbb{R}^d}  |  F _T(\bb{\xi}) |^2 d\bb{\xi}.$$

$\lim _{T \rightarrow \infty} |F _T(\bb{\xi})|^2$ can be interpreted as
the density of power per unit of frequency and is often referred to as
the *power spectral density* (PSD).

The same reasoning can be repeated for WSS random processes. While a
random process $\mathpzc{F}$ has infinite energy, it generally has
finite average power and one can define the PSD as 

$$\begin{aligned}
S _{\mathpzc{F}}(\bb{\xi}) =& \lim _{T \rightarrow \infty}  \frac{1}{T^\frac{d}{2}} \,  \mathbb{E} \left|\mathcal{F}\left( \mathrm{rect} _T \cdot \mathpzc{F} \right)   \right|^2  \\
=&
 \lim _{T \rightarrow \infty}    \frac{1}{T^d} \, \mathbb{E} \left( 
  \int _{[-\frac{T}{2},\frac{T}{2}]^d}  f(\bb{x}) e^{+2\pi \ii \, \bb{\xi}^\Tr \bb{x} } d\bb{x} 
\int _{[-\frac{T}{2},\frac{T}{2}]^d}  f(\bb{x}') e^{-2\pi \ii \, \bb{\xi}^\Tr \bb{x}' } d\bb{x}'
  \right) \\
  =& 
 \lim _{T \rightarrow \infty}    \frac{1}{T^d}   \int _{[-\frac{T}{2},\frac{T}{2}]^{2d}}  \mathbb{E} f(\bb{x})f(\bb{x}') \, e^{-2\pi \ii \, \bb{\xi}^\Tr (\bb{x}-\bb{x}') } d\bb{x}  d\bb{x}'.\end{aligned}$$

Changing the integration variable $\bb{x}'$ to
$\bb{y} = \bb{x}-\bb{x}'$, one obtains the following result known as the
*Wiener-Khinchin theorem*: 

$$\begin{aligned}
{S _{\mathpzc{F}} = \mathcal{F} R _{\mathpzc{F}}}\end{aligned}$$


This is a profound conclusion relating the PSD of a random signal to the
Fourier transform of its auto-correlation.

Since $R _{\mathpzc{F}}$ is an even function (i.e., invariant to space
mirroring), $S _{\mathpzc{F}}$ is real-valued. However, it is also easy
to show that $S _{\mathpzc{F}} \ge 0$. This stems from the fact that the
auto-correlation is positive semi-definite and, in case of wide-sense
stationarity, it is diagonalized by the Fourier transform; hence the its
spectrum is non-negative.

### Cross-spectral density

The result can be generalized to a pair of jointly WSS random signals
$\mathpzc{F}$ and $\mathpzc{G}$. We define the *cross-spectal density*
as

$$S _{\mathpzc{F}}(\bb{\xi}) = \lim _{T \rightarrow \infty}  \frac{1}{T^\frac{d}{2}} \, \mathbb{E} \left(  \mathcal{F}\left( \mathrm{rect} _T \cdot \mathpzc{F} \right)   )^\ast \mathcal{F}\left( \mathrm{rect} _T \cdot \mathpzc{G} \right)   
\right).$$ 

The Wiener-Kinchin theorem states 

$$\begin{aligned}
{S _{\mathpzc{F}\mathpzc{G}} = \mathcal{F} R _{\mathpzc{F}\mathpzc{G}}}\end{aligned}$$

where the cross-correlation is defined as
$R _{\mathpzc{F}\mathpzc{G}} (\bb{x}) = \mathbb{E} \mathpzc{F}(\bb{0}) \mathpzc{G}(\bb{x})$.

### LSI systems

Let $\mathpzc{F}$ be a WSS signal passing through a linear
shift-invariant system $\mathcal{H}$ with the impulse reponse $h$. We
define the output signal as $\mathpzc{G} = h \ast \mathpzc{F}$.
Straightforwardly, the mean function of $\mathpzc{G}$ is given by

$$\begin{aligned}
\mu _{\mathpzc{G}}(\bb{x}) =& \mathbb{E}  \mathpzc{G}(\bb{x}) =  \mathbb{E} \int _{\mathbb{R}^d} h(\bb{x}')  \mathpzc{F}(\bb{x} - \bb{x}') d\bb{x}' = 
 \int _{\mathbb{R}^d} h(\bb{x}')  \mu _\mathpzc{F}(\bb{x} - \bb{x}') d\bb{x}' \\
 =& \int _{\mathbb{R}^d} h(\bb{x}')   d\bb{x}'  \, \mu _\mathpzc{F}   = H(\bb{0}) \mu _\mathpzc{F},\end{aligned}$$

where $H(\bb{0})$ is usually called the *DC response* of $\mathcal{H}$.
Note that the mean function is constant.

The auto-correlation function of $\mathpzc{G}$ is given by

$$\begin{aligned}
R _{\mathpzc{G}}(\bb{x},\bb{x}-\bb{y}) =& \mathbb{E}  \mathpzc{G}(\bb{x}) \mathpzc{G}(\bb{x}+\bb{y}) \\
=&  \mathbb{E} \left( \int _{\mathbb{R}^d} h(\bb{x}')  \mathpzc{F}(\bb{x} - \bb{x}') d\bb{x}'   \int _{\mathbb{R}^d} h(\bb{x}'')  \mathpzc{F}(\bb{x}-\bb{y} - \bb{x}'') d\bb{x}''  \right) \\
=&   \int _{\mathbb{R}^d} h(\bb{x}'')   \left( \int _{\mathbb{R}^d} h(\bb{x}') \,  \mathbb{E}   \mathpzc{F}(\bb{x} - \bb{x}') \mathpzc{F}(\bb{x}-\bb{y} - \bb{x}'') d\bb{x}' \right) d\bb{x}'' \\
=&   \int _{\mathbb{R}^d} h(\bb{x}'')   \left( \int _{\mathbb{R}^d} h(\bb{x}') \,  R _\mathpzc{F} (\bb{y} + \bb{x}'' - \bb{x}') d\bb{x}' \right) d\bb{x}'' \\
=&   \int _{\mathbb{R}^d} h(\bb{x}'')   (h \ast R _\mathpzc{F}) (\bb{y} + \bb{x}'')  d\bb{x}'' 
= \int _{\mathbb{R}^d} \overline{h}(\bb{y}'')   (h \ast R _\mathpzc{F}) (\bb{y} - \bb{y}'')  d\bb{y}''
\\
=& (h \ast R _\mathpzc{F} \ast \overline{h})(\bb{y}),\end{aligned}$$

where $ \overline{h}(\bb{x}) = h(-\bb{x})$. Note that the
auto-correlation function is shift-invariant (i.e., it does not depend
on $\bb{x}$). These two results imply that $\mathpzc{Y}$ is WSS (not
surprising: a signal with shift-invariant moments passing through a
shift-invariant system leads to a signal with shift-invariant moments).

The input-output cross-correlation is given by 

$$\begin{aligned}
R _{\mathcal{F}\mathpzc{G}}(\bb{x},\bb{x}-\bb{y}) =& \mathbb{E}  \mathpzc{F}(\bb{x}) \mathpzc{G}(\bb{x}-\bb{y}) \\
=&  \mathbb{E} \left( \mathpzc{F}(\bb{x})   \int _{\mathbb{R}^d} h(\bb{x}')  \mathpzc{F}(\bb{x}-\bb{y} - \bb{x}') d\bb{x}' \right) \\
=&   \int _{\mathbb{R}^d} h(\bb{x}') \,  \mathbb{E}   \mathpzc{F}(\bb{x}) \mathpzc{F}(\bb{x}-\bb{y} - \bb{x}') d\bb{x}'  \\
=&   \int _{\mathbb{R}^d} h(\bb{x}') \,  R _\mathpzc{F}(\bb{y} + \bb{x}') d\bb{x}'  = (\bar{h} \ast R _\mathpzc{F})(\bb{y}).\end{aligned}$$


Again note that the function is shift-invariant, implying that
$\mathpzc{F}$ and $\mathpzc{G}$ are jointly WSS. Clearly, since
$R _{\mathcal{G}\mathpzc{F}} = \overline{R} _{\mathcal{F}\mathpzc{G}} = {h} \ast \overline{R} _\mathpzc{F}  = {h} \ast {R} _\mathpzc{F}$.

Translating the latter expressions to the frequency domain, we obtain
the following relations 

$$\begin{aligned}
R _{\mathpzc{G}} = h \ast R _{\mathpzc{F}} \ast \overline{h} & \,\, \fff \,\, & S _{\mathpzc{G}} = |H|^2 \cdot S _{\mathpzc{F}} \\
R _{\mathpzc{F} \mathpzc{G}} = \overline{h} \ast R _{\mathpzc{F}}  & \,\, \fff \,\,& S _{\mathpzc{F}\mathpzc{G}} = H^\ast \cdot S _{\mathpzc{F}} \\
R _{\mathpzc{G} \mathpzc{H}} = {h} \ast R _{\mathpzc{F}}  & \,\, \fff \,\,& S _{\mathpzc{G}\mathpzc{F}} = H \cdot S _{\mathpzc{F}}
\end{aligned}
$$ 

Compare this to linear transformations of random vectors.

### White and colored noise

A random signal $\mathpzc{N}$ with constant PSD,
$S _{\mathpzc{N}}(\bb{\xi}) = \sigma _\mathpzc{N}^2$ is usually called
*white noise*, by analogy with white light that has approximately flat
spectrum[^3]. Its auto-correlation is given by
$R _{\mathpzc{N}}(\bb{x}) = \sigma _\mathpzc{N}^2 \, \delta$. When white
noise passes through an LSI system $\mathcal{H}$, the spectrum at the
output,
$S _{\mathcal{H}\mathpzc{N}}(\bb{\xi}) = |H(\bb{\xi})|^2 \sigma _\mathpzc{N}^2$,
is shaped by the power response $|H(\bb{\xi})|^2$ of the system. This
phenomenon is called as *coloring*. The auto-correlation function of
colored noise is given by
$R _{\mathcal{H}\mathpzc{N}} = \sigma _\mathpzc{N}^2 \, h \ast  \delta \ast \overline{h} = \sigma _\mathpzc{N}^2 \, h \ast \overline{h}$.

[^1]: To be completely rigorous, the proper way to define the PDF is by
    first equipping the image of the map $\mathpzc{X}$ with the Lebesgue
    measure $\lambda$ that assigns to every interval $[a,b]$ its length
    $b-a$. Then, we invoke the Radon-Nikodym theorem saying that if
    $\mathpzc{X}$ is absolutely continuous w.r.t. $\lambda$, there
    exists a measurable function $f : \mathbb{R} \rightarrow [0,\infty)$
    such that for every measurable $A \subset \mathbb{R}$,
    $\displaystyle{(\mathpzc{X} _\ast P)(A) =P(\mathpzc{X}^{-1}(A)) = \int _A f d\lambda}$.
    $f$ is called the *Radon-Nikodym derivative* and denoted by
    $\displaystyle{f = \frac{d(\mathpzc{X} _\ast P)}{d\lambda}} $. It is
    exactly our PDF.

[^2]: In fact, $r _{\mathpzc{X}\mathpzc{Y}}$ can be viewed as an inner
    product on the space of random variables. This creates a geometry
    isomorphic to the standard Euclidean metric in $\mathbb{R}^n$. Using
    this construction, the Cauchy-Schwarz inequality immediately
    follows:
    $| r _{\mathpzc{X}\mathpzc{Y}} | \le \sigma _\mathpzc{X} \sigma _\mathpzc{Y}$.

[^3]: This appears to be false when color perception is examined more
    closely!
