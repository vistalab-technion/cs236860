---
title: "Lecture 2: Sampling and Interpolation"
excerpt: Dirac's delta, Sampling Theorem, Interpolation
author: alex
copyright: alex
---

## Notation

In our course, we will deal almost exclusively with real functions. A
scalar function will be denoted as $f : \RR^d \rightarrow \RR$,
$f({\bm{\mathrm{x}}})$, or simply $f$, with $d$ denoting the domain
dimension. $d$-dimensional multi-indices will be denoted by
${\bm{\mathrm{n}}} = (n_1,\dots,n_d) \in {\bm{\mathrm{Z}}}^d$. An
operator acting on a function will be denoted by $\mathcal{H}(f)$ or
simply $\mathcal{H}f$.

## Dirac’s delta

In what follows, we define a little bit of a technical machinery that
will have crucial importance in the rigorous treatment of sampling.

### Adjoint operators

We start with a little digression. Let $\mathbb{U}$ and $\mathbb{V}$ be
two spaces equipped with inner products, and let
$\mathcal{A} : \mathbb{U} \rightarrow \mathbb{V}$ and
$\mathcal{B} : \mathbb{V} \rightarrow \mathbb{U}$ be two operators. If
for every $u \in \mathbb{U}$ and $v \in \mathbb{V}$ it holds that
$\langle \mathcal{A}u, v \rangle_{\mathbb{V}} = \langle u, \mathcal{B} v \rangle_{\mathbb{U}}$,
we will say that the operator $\mathcal{B}$ is the *adjoint* of
$\mathcal{A}$ and denote it by $\mathcal{B} = \mathcal{A}^\ast$. Adjoint
is the generalization of the notion of transpose (Hermitian transpose in
the complex case) of a matrix. An operator satisftying
$\mathcal{A}^\ast  = \mathcal{A}$ is said to be *self-adjoint* (think of
a symmetric matrix). Warning: it is tempting to confuse adjoint with
inverse, as both are operators from $\mathbb{V}$ to $\mathbb{U}$. These
are very different notions, coinciding only in the case of unitary
operators.

For example, the adjoint of the translation operator
$\tau_{\bm{\mathrm{p}}}$ can be immediately derived from

$$
\begin{aligned}
\langle \tau_{\bm{\mathrm{p}}}f, g \rangle = \int_{\RR^d} f(\bm{\mathrm{x}} - \bm{\mathrm{p}}) g^\ast (\bm{\mathrm{x}}) d\bm{\mathrm{x}} = \int_{\RR^d} f(\bm{\mathrm{y}}) g^\ast (\bm{\mathrm{y}}-(-\bm{\mathrm{p}}) ) d\bm{\mathrm{y}}  = \langle f, \tau_{-\bm{\mathrm{p}}}g \rangle,
\end{aligned}
$$

from where $\tau_{\bm{\mathrm{p}}}^\ast = \tau_{-\bm{\mathrm{p}}}$
(verify this by thinking of $\tau_{\bm{\mathrm{p}}}$ as of a Toeplitz
matrix and taking its transpose). As a more general example, consider
operator convolving the input function with a given real-valued function
$h$. Then, 

$$
\begin{aligned}
\langle h \ast f, g \rangle &= \int_{\RR^d} \left(  \int_{\RR^d} h(\bm{\mathrm{x}}-\bm{\mathrm{x}}') f(\bm{\mathrm{x}}') d\bm{\mathrm{x}}' \right) g^\ast (\bm{\mathrm{x}}) d\bm{\mathrm{x}}  \\
&= \int_{\RR^d}f(\bm{\mathrm{x}}') \left(  \int_{\RR^d} h(\bm{\mathrm{x}}- \bm{\mathrm{x}}') g(\bm{\mathrm{x}}) d\bm{\mathrm{x}} \right)^\ast d\bm{\mathrm{x}}' \\
&= \int_{\RR^d}f(\bm{\mathrm{x}}') \left(  \int_{\RR^d} \bar{h}(\bm{\mathrm{x}}'-\bm{\mathrm{x}}) g(\bm{\mathrm{x}}) d\bm{\mathrm{x}} \right)^\ast d\bm{\mathrm{x}}' = \langle f, \bar{h} \ast g \rangle,
\end{aligned}
$$

from where the adjoint of $h \ast f$ is $\bar{h} \ast f$, with
$\bar{h}({\bm{\mathrm{x}}}) = h(-{\bm{\mathrm{x}}})$.

Yet another example: it is straightforward to see that the point-wise
product with a given real-valued function $h$ is self-adjoint, since

$$
\begin{aligned}
\langle h \cdot f, g \rangle &= \int_{\RR^d} h({\bm{\mathrm{x}}}) f({\bm{\mathrm{x}}}) g^\ast ({\bm{\mathrm{x}}}) d{\bm{\mathrm{x}}}\\
 &= \int_{\RR^d}  f({\bm{\mathrm{x}}}) h^\ast ({\bm{\mathrm{x}}}) g^\ast ({\bm{\mathrm{x}}}) d{\bm{\mathrm{x}}}\\
 &= \langle f, h \cdot g \rangle.
\end{aligned}
$$

Finally, let us derive the adjoint of the Fourier transform:

$$
\begin{aligned}
\langle \mathcal{F} f, G \rangle &= \int_{\RR^d}  \left(  \int_{\RR^d} f({\bm{\mathrm{x}}}) e^{-2\pi \ii \, {\bm{\mathrm{\xi}}}^\Tr {\bm{\mathrm{x}}}}  d{\bm{\mathrm{x}}} \right) G^\ast ({\bm{\mathrm{\xi}}}) d{\bm{\mathrm{\xi}}}\\ 
&= \int_{\RR^d}  f({\bm{\mathrm{x}}}) \left(  \int_{\RR^d} G ({\bm{\mathrm{\xi}}})  e^{2\pi \ii \, {\bm{\mathrm{\xi}}}^\Tr {\bm{\mathrm{x}}}}  d{\bm{\mathrm{\xi}}}  \right)^\ast d{\bm{\mathrm{x}}} \\
&= \langle  f, \mathcal{F}^{-1} G \rangle.
\end{aligned}
$$

We observe that $\mathcal{F}^\ast = \mathcal{F}^{-1}$, implying that $\mathcal{F}$
is a unitary transformation (an isometry). From this geometric
properties, *Plancherel’s identity* follows:

$$
\langle f, g \rangle = \langle F, G \rangle,
$$ 

with
$F = \mathcal{F} f$ and $G = \mathcal{F} g$. In the particular case of
$g=f$, we obtain the celebrate *Parseval’s identity*,

$$
\|f\|^2 = \langle f, f \rangle = \langle F, F \rangle = \| F \|^2,
$$

which is precisely what *isometry* means.

We will see more adjoint operators in the sequel.

### Distributions

Let $f : \RR^d \rightarrow \RR$ be given function. A funny way to
represent it is to consider the action of the inner product with $f$ on
a real-valued test function[^1] $\psi$:

$$
F(\psi) = \langle f, \psi \rangle = \int_{\RR^d} f({\bm{\mathrm{x}}}) \psi({\bm{\mathrm{x}}}) d{\bm{\mathrm{x}}}.
$$

Essentially, we defined a linear functional assigning each test function
a real scalar. Let us now think of a probability measure $\mu$ on
$\RR^d$. A way to represent the measure is to consider its moments of
the form $$M(\psi) = \mathbb{E}_\mu \psi = \int \psi d\mu$$ (we can
think of the previous example as a realization of the latter one if we
think of $f$ as the probability density function embodying $\mu$).
Again, a linear functional is used to represent a probability measure.



A general linear functional assigning real scalars to tests functions is
called a *distribution*. Because of our first example, we will denote
the action of such a functional on a test function using the inner
product notation, $$L(\psi) = \langle L, \psi \rangle.$$ Beware: this is
just a convenient notation. Of course, $L$ is not a function and the
above inner product makes no sense. We simply *define* it to assume the
value of $L(\psi)$. The notation is convenient because it allows to
define in a consistent way what it means to apply an operator to a
distribution. 


Let us return to our first example of the distribution
$F(\psi) = \langle f, \psi \rangle$ representing a function $f$, and let
us now be given an operator $\mathcal{A}$. What is the meaning of
“$\mathcal{A}F$”? Rigorously, $F$ is not a function and the latter
expression makes no sense; yet, can we *attribute* it a meaning? In our
particular case, we can interpret $\mathcal{A}F$ as a new distribution

$$
\mathcal{A}F)(\psi) = \langle \mathcal{A} f, \psi \rangle = \langle  f, \mathcal{A}^\ast \psi \rangle.
$$

We will now carry this trick with the adjoint to the general case: given
a distribution $L$, we *define*

$$
(\mathcal{A} L)(\psi) = \langle \mathcal{A} L, \psi \rangle = \langle  L, \mathcal{A}^\ast \psi \rangle = L(\mathcal{A}^\ast \psi).
$$


Note that the left hand side is just convenient notation that rigorously
speaking has no mathematical sense. The right hand side, on the contrary
is well-defined: no operator is acting on $L$; instead, the adjoint
operator acts on the test function. This essentially gives rigorous
sense and meaning to “$\mathcal{A} L$”. With these understandings in
mind, let us now introduce a very important distribution and study its
properties.

### Delta distribution

The distribution $\delta(\psi) = \psi(0)$ is called the *Dirac delta*.
Using the inner product notation, we can write

$$
\delta(\psi) = \langle \delta, \psi \rangle = \int_{\RR^d} \delta({\bm{\mathrm{x}}}) \psi(x) d{\bm{\mathrm{x}}}.
$$

While the second expression is just a syntactic sugar to $\delta(\psi)$,
the third one is simply nonsensensical: $\delta$ is not a function and
cannot be integrated; furthermore, there exists no such a function the
integration with which would yield the value of $\psi$ at
${\bm{\mathrm{x}}} = 0$.[^2] ! Yet, engineering books are full with such
a notation, referring to $\delta$ as to a “generalized function”.

Using our definition of
$(\mathcal{A} \delta)(\psi) = \delta(\mathcal{A}^\ast \psi)$, we can
readily list a few important properties of the delta (it only takes to
derive the adjoint operator).

### Translation

$$
(\tau_{\bm{\mathrm{p}}} \delta)(\psi) = \delta( \tau_{\bm{\mathrm{p}}}^\ast \psi) = \delta( \tau_{-\bm{\mathrm{p}}} \psi) = \left. \psi(\bm{\mathrm{x}} + \bm{\mathrm{p}}) \right| _{\bm{\mathrm{x}}=0} = \psi(\bm{\mathrm{p}}).
$$

### Pointwise product

Given a function $h : \RR^d \rightarrow \RR$ and using the fact that
pointwise product with a real-valued function is self-adjoint, we obtain

$$
(h \cdot \delta)(\psi) = \delta( h \cdot \psi) = f(0) \psi(0) = h(0)\, \delta(\psi).
$$

### Convolution

Given a function $h : \RR^d \rightarrow \RR$ and using the fact that the
adjoint of convolution with $h$ is the convolution with
$\bar{h}({\bm{\mathrm{x}}}) = h(-{\bm{\mathrm{x}}})$ yields

$$
(h \ast \delta)(\psi) = \delta( \bar{h} \ast \psi) =   (\bar{h} \ast \psi)(0) = \int_{\RR^d} h(\bm{\mathrm{x}}-0) \psi (\bm{\mathrm{x}}) d\bm{\mathrm{x}} = \langle h, \psi \rangle.
$$


Note that the right hand side is a distribution representing the
function $h$. With some abuse of notation, we can say that the
convolution of $\delta$ with a function $h$ ceases to be distribution
and becomes the function $h$ itself. Delta is actually the only
distribution acting as an identity element, a fact allowing to define a
group of functions with convolution serving as the group operation.

### Fourier transform

Using the unitarity of the Fourier transform, 

$$
\begin{aligned}
(\mathcal{F} \delta)(\Psi) &=& ( \delta)(\mathcal{F}^\ast \Psi) = (\mathcal{F}^{-1} \Psi)(0) = \int_{\RR^d} \Psi({\bm{\mathrm{\xi}}}) e^{ 2\pi \ii \, {\bm{\mathrm{\xi}}}^\Tr {\bm{\mathrm{0}}} } d{\bm{\mathrm{\xi}}}
=  \int_{\RR^d} \Psi({\bm{\mathrm{\xi}}}) \mathbbl{1}({\bm{\mathrm{\xi}}}) d{\bm{\mathrm{\xi}}} = \langle \Psi, \mathbbl{1} \rangle,
\end{aligned}
$$

where $\mathbbl{1}({\bm{\mathrm{\xi}}}) = 1$ is a constant function.
Again, the last expression is a distribution representing the constant
function $\mathbbl{1}$, so we can write (with some abuse of notation):
$ \delta  \fff \mathbbl{1}$. Combining this result with the
translation/modulation property leads to 

$$
\begin{aligned}
{\tau_{\bm{\mathrm{p}}} \delta \fff \phi_{\bm{\mathrm{p}}} }
\end{aligned}
$$

### Stretching

In order to see the action of the stretrching operator
$\mathcal{S}_{\bm{\mathrm{A}}} : f(\bm{\mathrm{x}}) \mapsto f(\bm{\mathrm{A}}\bm{\mathrm{x}})$ on the delta, let us first derive its adjoint: 

$$
\begin{aligned}
\langle \mathcal{S}_{\bm{\mathrm{A}}} f, g \rangle &= \int_{\RR^d} f(\bm{\mathrm{A}} \bm{\mathrm{x}}) g(\bm{\mathrm{x}}) d\bm{\mathrm{x}}\\
 &= \int_{\RR^d} f(\bm{\mathrm{y}}) g(\bm{\mathrm{A}}^{-1} \bm{\mathrm{y}}) \frac{d\bm{\mathrm{y}}}{| \det \bm{\mathrm{A}} |} \\
 &= \langle f , \frac{ \mathcal{S}_{\bm{\mathrm{A}}}^{-1} }{ | \det {\bm{\mathrm{A}}} | } g \rangle,
\end{aligned}
$$

from where
$ {\mathcal{S} _\bm{\mathrm{A}}}^\ast = \frac{ \mathcal{S} _{\bm{\mathrm{A}}}^{-1} }{ | \det {\bm{\mathrm{A}}} | } $.

Invoking the definion of an operator acting on a distribution,

$$
( \mathcal{S}_{\bm{\mathrm{A}}} \delta)(\psi) = \delta( \mathcal{S}_{\bm{\mathrm{A}}}^\ast \psi) = \left. \frac{\psi(\bm{\mathrm{A}}^{-1} \bm{\mathrm{x}})}{| \det \bm{\mathrm{A}} |} \right|_{\bm{\mathrm{x}} =0} =  \frac{\psi(0)}{| \det \bm{\mathrm{A}} |}   =  \frac{ \delta(\psi) }{| \det \bm{\mathrm{A}} |}.
$$	

Informally, thinking of $\delta$ as of a “function”, we can write
$\delta({\bm{\mathrm{A}}}{\bm{\mathrm{x}}})  = \frac{ \delta({\bm{\mathrm{x}}}) }{| \det {\bm{\mathrm{A}}} |}$.
A straightforward corollary is that $\delta$ is invariant to rotation.

### Derivatives

In order to define delta’s derivatives, let us first examine the adjoint
of the derivative operator. For the sake of clairity, let us first
consider the partial derivative with respect to the first coordinate,
$\partial_1 f = \frac{\partial}{\partial x_1} f$. Using integration by
parts yields 

$$
\begin{aligned}
\langle \partial_1 f, g \rangle &= \int_{\RR^d} \partial_1 f({\bm{\mathrm{x}}}) g({\bm{\mathrm{x}}}) d{\bm{\mathrm{x}}} \\
&= \left. \int_{\RR^{d-1}}  f({\bm{\mathrm{x}}}) g({\bm{\mathrm{x}}}) dx_2 \cdots dx_d \right|_{x_1= -\infty}^{x_1 = \infty} -   \int_{\RR^d}  f({\bm{\mathrm{x}}})\partial_1 g({\bm{\mathrm{x}}}) d{\bm{\mathrm{x}}}.
\end{aligned}
$$


By asserting that the functions $f$ and $g$ vanish at infinity, we
obtain 

$$\begin{aligned}
\langle \partial_1 f, g \rangle = \langle  f, - \partial_1 g \rangle.
\end{aligned}$$

In the same manner, for any partial derivative operator of order $n$,

$$
\partial_{\bm{\mathrm{n}}} f = \frac{\partial^{n_1}}{\partial x_1^{n_1}} \cdots \frac{\partial^{n_d}}{\partial x_d^{n_d}} f,
$$

parametrized by the multi-index ${\bm{\mathrm{n}}}$, it is
straightforward to show that
$\partial_{\bm{\mathrm{n}}}^\ast = (-1)^{|\bm{\mathrm{n}}|} \partial_{\bm{\mathrm{n}}}$,
where $|\bm{\mathrm{n}}| = n_1 + \cdots + n_d$. With this result, we
can now write

$$
( \partial_{\bm{\mathrm{n}}}  \delta)(\psi) = \delta( \partial_{\bm{\mathrm{n}}}^\ast \psi) = (-1)^{|\bm{\mathrm{n}}|} (\partial_{\bm{\mathrm{n}}} \psi )(0).
$$

## Dirac’s bed


### Dirac’s bed

Equipped with the delta and its basic properties, we are ready to define
a new distribution that we will call a *Dirac’s bed* or a *bed of
impulses* (on the integer lattice; we will generalize this definition in
the sequel):

$$
\bed = \sum_{\bm{\mathrm{p}} \in \mathbb{Z}^d} \tau_{\bm{\mathrm{p}}} \delta.
$$

The importance of this distribution comes from the following two
properties:

### Pointwise product

Given a real-valued function $f$, 

$$
\begin{aligned}
(f \cdot \bed)(\psi) &=  \sum_{\bm{\mathrm{p}} \in \mathbb{Z}^d} \tau_{\bm{\mathrm{p}}} (f \cdot \delta)(\psi)
= \sum_{\bm{\mathrm{p}} \in \mathbb{Z}^d}  (\tau_{\bm{\mathrm{p}}} \delta)(f \cdots \psi) = 
\sum_{\bm{\mathrm{p}} \in \mathbb{Z}^d}   \delta(\tau_{-\bm{\mathrm{p}}} (f \cdot \psi)) = \sum_{\bm{\mathrm{p}} \in \mathbb{Z}^d}   f({\bm{\mathrm{p}}}) \psi(\bm{\mathrm{p}}) \nonumber\\
&=  \sum_{\bm{\mathrm{p}} \in \mathbb{Z}^d}   f(p) (\tau_{\bm{\mathrm{p}}} \delta)(\psi).
\end{aligned}
$$

Therefore, we can write 

$$
\begin{aligned}
 f \cdot \bed = \sum_{\bm{\mathrm{p}} \in \mathbb{Z}^d}   f(p) \, \tau_{\bm{\mathrm{p}}} \delta  
\end{aligned}
$$

Note that the resulting distribution is still a bed of impulses, but now
each impulse placed at point ${\bm{\mathrm{p}}}$ of the integer lattice
is scaled by the value of the function $f$ at that point. This is a
convenient way to represent the action of *sampling*.

### Convolution

As in the previous case, 

$$
\begin{aligned}
(f \ast \bed)(\psi) &=  \sum_{\bm{\mathrm{p}} \in \mathbb{Z}^d} \tau_{\bm{\mathrm{p}}} (f \ast \delta)(\psi)
= \sum_{\bm{\mathrm{p}} \in \mathbb{Z}^d}  (\tau_{\bm{\mathrm{p}}} \delta)(\bar{f} \ast \psi) = 
\sum_{\bm{\mathrm{p}} \in \mathbb{Z}^d}   (\bar{f} \ast \psi)(\bm{\mathrm{p}}) \\
&= \sum_{\bm{\mathrm{p}} \in \mathbb{Z}^d} \int_{\RR^d} f(\bm{\mathrm{x}}-\bm{\mathrm{p}}) g(\bm{\mathrm{x}}) d\bm{\mathrm{x}}\\
&= \int_{\RR^d} \left( \sum_{\bm{\mathrm{p}} \in \mathbb{Z}^d}  f(\bm{\mathrm{x}}-\bm{\mathrm{p}}) \right) g(\bm{\mathrm{x}}) d\bm{\mathrm{x}} \\
&= \langle \sum_{\bm{\mathrm{p}} \in \mathbb{Z}^d}  f(\bm{\mathrm{x}}-\bm{\mathrm{p}}), g \rangle.
\end{aligned}
$$

Note that the resulting distribution is no more a bed of impulses but
rather the periodic function 

$$
\begin{aligned}
f \ast \bed =\sum_{\bm{\mathrm{p}} \in \mathbb{Z}^d}  f(\bm{\mathrm{x}}-\bm{\mathrm{p}}) 
\end{aligned}
$$


In other words, convolution with a bed of impulses is akin to
*periodization* of the function.

### Poisson summation identity

Let $\psi$ be a test function and let us define it periodized version

$$
\tilde{\psi}(\bm{\mathrm{x}}) = (\bed \ast \psi)(\bm{\mathrm{x}}) = \sum_{\bm{\mathrm{p}} \in \mathbb{Z}^d}  \psi(\bm{\mathrm{x}}-\bm{\mathrm{p}}).
$$

As a periodic function, it can be represented as a ($d$-dimensional)
Fourier series

$$
\tilde{\psi}(\bm{\mathrm{x}}) = \sum_{\bm{\mathrm{p}} \in \mathbb{Z}^d} c_{\bm{\mathrm{p}}}  e^{2\pi \ii \, \bm{\mathrm{p}}^\Tr \bm{\mathrm{x}}}
$$

with 

$$
\begin{aligned}
c_{\bm{\mathrm{p}}} &= \int_{[0,1]^d}  \tilde{\psi}(\bm{\mathrm{x}})  e^{-2\pi \ii \, \bm{\mathrm{p}}^\Tr \bm{\mathrm{x}}} d\bm{\mathrm{x}} =   \int_{[0,1]^d} \sum_{\bm{\mathrm{q}} \in \mathbb{Z}^d}  \psi(\bm{\mathrm{x}}-\bm{\mathrm{q}})  e^{-2\pi \ii \, \bm{\mathrm{p}}^\Tr \bm{\mathrm{x}}} d\bm{\mathrm{x}} \\
&= \sum_{\bm{\mathrm{q}} \in \mathbb{Z}^d}   \int_{[0,1]^d} \psi(\bm{\mathrm{y}}) e^{-2\pi \ii \, \bm{\mathrm{p}}^\Tr (\bm{\mathrm{y}} + \bm{\mathrm{q}})} d\bm{\mathrm{y}}  \\
&= \sum_{\bm{\mathrm{q}} \in \mathbb{Z}^d}   \int_{[0,1]^d} \psi(\bm{\mathrm{y}}) e^{-2\pi \ii \, \bm{\mathrm{p}}^\Tr \bm{\mathrm{y}} } d \bm{\mathrm{y}}   \, e^{2\pi \ii \, \bm{\mathrm{p}}^\Tr \bm{\mathrm{q}}} \\
&= \int_{\RR^d} \psi(\bm{\mathrm{y}}) e^{-2\pi \ii \, \bm{\mathrm{p}}^\Tr \bm{\mathrm{y}} } d\bm{\mathrm{y}} = (\mathcal{F}\psi)(\bm{\mathrm{p}}).
\end{aligned}
$$

Therefore,

$$
\tilde{\psi}(\bm{\mathrm{x}}) = \sum_{\bm{\mathrm{p}} \in \mathbb{Z}^d}  (\mathcal{F}\psi)(\bm{\mathrm{p}})  e^{2\pi \ii \, \bm{\mathrm{p}}^\Tr \bm{\mathrm{x}}}.
$$

Substituting $\bm{\mathrm{x}} = \bm{\mathrm{0}}$, we obtain

$$
\tilde{\psi}(\bm{\mathrm{0}}) =  \sum_{\bm{\mathrm{p}} \in \mathbb{Z}^d}  \psi(-\bm{\mathrm{p}}) = \sum_{\bm{\mathrm{p}} \in \mathbb{Z}^d}  \psi(\bm{\mathrm{p}})
$$

as well as

$$
\tilde{\psi}(\bm{\mathrm{0}}) =  \sum_{\bm{\mathrm{p}} \in \mathbb{Z}^d}  (\mathcal{F}\psi)(\bm{\mathrm{p}})  e^{2\pi \ii \, \bm{\mathrm{p}}^\Tr \bm{\mathrm{0}}} = \sum_{\bm{\mathrm{p}} \in \mathbb{Z}^d}  (\mathcal{F}\psi)(\bm{\mathrm{p}}).
$$

This leads to the following amazing result known as the *Poisson
summation formula*:

$$
\sum_{\bm{\mathrm{p}} \in \mathbb{Z}^d}  (\mathcal{F}\psi)(\bm{\mathrm{p}}) = \sum_{\bm{\mathrm{p}} \in \mathbb{Z}^d}  \psi(\bm{\mathrm{p}}).
$$

The following corollary follows directly: 

$$
\begin{aligned}
\langle \mathcal{F} \bed, \Psi \rangle &= \langle  \bed,  \mathcal{F}^{-1} \Psi \rangle = \langle  \bed,  \psi \rangle  
= \sum_{\bm{\mathrm{p}} \in \mathbb{Z}^d}  \psi(\bm{\mathrm{p}}) = \sum_{\bm{\mathrm{p}} \in \mathbb{Z}^d}  \Psi({\bm{\mathrm{p}}}) = \langle \bed, \Psi \rangle
\end{aligned}
$$

where we denote as usual $\Psi = \mathcal{F} \psi$. In other words,

$$
\begin{aligned}
  \mathcal{F} \bed = \bed 
\end{aligned}
$$

Combining this result with the convolution/product duality, we obtain
the following duality between sampling and periodization:

$$
\begin{aligned}
{ \bed \ast f \fff \bed \cdot F }
\end{aligned}
$$

In other words, sampling the signal in the space domain result in the
periodization of the Fourier transform, and vice versa, periodization of
the signal in the space domain results in a sampled Fourier transform.
The latter is a mere justification of the fact that periodic functions
can be represented as Fourier series that have a discrete set of spatial
frequencies. The former leads to the famous sampling theorem.

## Sampling theorem


Sampling a function on the integer lattice is descrbied by the following
distribution

$$
\bed \cdot f = \sum_{\bm{\mathrm{p}} \in \mathbb{Z}^d} f(\bm{\mathrm{p}}) \tau_{\bm{\mathrm{p}}} \delta.
$$

Our previous result tells us that the Fourier transform of this
distribution is the periodization of the Fourier transform
$F({\bm{\mathrm{\xi}}}) = (\mathcal{F}f)({\bm{\mathrm{\xi}}})$,

$$
\tilde{F}(\bm{\mathrm{\xi}}) =(\bed \ast F)(\bm{\mathrm{\xi}}) = \sum_{\bm{\mathrm{p}} \in \mathbb{Z}^d}  F(\bm{\mathrm{\xi}} - \bm{\mathrm{p}}).
$$

It it happens that
$\mathrm{supp}(F) = \{  {\bm{\mathrm{\xi}}} : F({\bm{\mathrm{\xi}}}) \ne 0 \} \not\subset \left[-\frac{1}{2},\frac{1}{2}\right]^d$,
in the above summation, the shifted replicas of the Fourier transforms
of $f$ will overlap forever losing the information contained in them and
producing originally inexistent content at lower frequencies. This
phenomenon is called *aliasing*.

On the other hand, if $\mathrm{supp}(F) \subset \left[-\frac{1}{2},\frac{1}{2}\right]^d$, the
shifted replicas will not overlap and no information loss will occur.

The original Fourier transform of $f$ can be recovered from $\tilde{F}$
by multiplying the latter with the box function taking the value of $1$
on $\left[-\frac{1}{2},\frac{1}{2}\right]^d$ and $0$ elsewhere:

$$
F = \tilde{F} \cdot \mathrm{box}
$$. 

We have seen already the following duality $$\mathrm{box} \fff \mathrm{sinc},$$ from which we obtain

$$
f(\bm{\mathrm{x}}) = (\bed \cdot f) \ast \mathrm{sinc} (\bm{\mathrm{x}}) =  \sum_{\bm{\mathrm{p}} \in \mathbb{Z}^d} f(\bm{\mathrm{p}}) \tau_{\bm{\mathrm{p}}} \delta \ast \mathrm{sinc}(\bm{\mathrm{x}}) = \sum_{\bm{\mathrm{p}} \in \mathbb{Z}^d} f(\bm{\mathrm{p}}) \, \mathrm{sinc}(\bm{\mathrm{x}} - \bm{\mathrm{p}}).
$$

This result tells us that if the function $f$ has
$\mathrm{supp}(F) \subset \left[-\frac{1}{2},\frac{1}{2}\right]^d$ (such
functions are called *band limited*), we can perfectly reconstruct it
from its sampled version by means of the above series (which is
sometimes called the *cardinal series*), and is known as the
(Kotelnikov-Whittaker-Nyquist-Hartley-Shannon) *sampling theorem*.

In what follows, we extend this result to sampling on a general lattice.

### Lattices

The sampling and periodization we dealt with so far was defined on the
so-called *integer lattice*, $\mathbb{Z}^d$. A more general lattice is
defined by means of a regular (=non-singular) linear transformation of
$\mathbb{Z}^d$:

$$
\mathcal{L} = {\bm{\mathrm{A}}}\mathbb{Z}^d = \{  z_1 {\bm{\mathrm{a}}}_1 + \cdots + z_d {\bm{\mathrm{a}}}_d : {\bm{\mathrm{z}}} \in \mathbb{Z}^d \}.
$$

To perform sampling on $\mathcal{L} = {\bm{\mathrm{A}}}\mathbb{Z}^d$, we
define

$$
\bed_{\mathcal{L}} = \sum_{\bm{\mathrm{p}} \in \mathcal{L}} \tau_{\bm{\mathrm{p}}} \delta = \sum_{\bm{\mathrm{q}} \in \mathbb{Z}^d}  \tau_{\bm{\mathrm{A}}\bm{\mathrm{q}}} \delta.
$$

Since

$$ 
(\tau_{\bm{\mathrm{A}}\bm{\mathrm{q}}} f)(\bm{\mathrm{x}}) = f(\bm{\mathrm{x}}-\bm{\mathrm{A}}\bm{\mathrm{q}}) = f(\bm{\mathrm{A}}(\bm{\mathrm{A}}^{-1}(\bm{\mathrm{x}}-\bm{\mathrm{A}}\bm{\mathrm{q}}))) = ((\mathcal{S}_{\bm{\mathrm{A}}^{-1}} \tau_{\bm{\mathrm{q}}}  \mathcal{S}_{\bm{\mathrm{A}}}) f)(\bm{\mathrm{x}}) 
$$,

we have

$$
\bed_{\mathcal{L}} = \sum_{\bm{\mathrm{q}} \in \mathbb{Z}^d}  \mathcal{S}_{\bm{\mathrm{A}}^{-1}} \tau_{\bm{\mathrm{q}}} \mathcal{S}_{\bm{\mathrm{A}}}  \delta 
= \frac{\mathcal{S}_{\bm{\mathrm{A}}^{-1}}}{ |\det \bm{\mathrm{A}} |}  \sum_{\bm{\mathrm{q}} \in \mathbb{Z}^d}   \tau_{\bm{\mathrm{q}}} \delta =  \frac{\mathcal{S}_{\bm{\mathrm{A}}^{-1}}}{ \mathrm{vol}\, \mathcal{L} } \bed,
$$

where by $\mathrm{vol}\, \mathcal{L}$ we denote the volume of the
structural element of the lattice.

Using the stretching property of the Fourier transform,

$$
\mathcal{S}_{\bm{\mathrm{A}}} \fff \frac{\mathcal{S}_{\bm{\mathrm{A}}^{-\Tr}}}{ | \det \bm{\mathrm{A}} | },
$$

we have

$$
\begin{aligned}
\mathcal{F} \bed_{\mathcal{L}} &= \frac{1}{ | \det \bm{\mathrm{A}} | } \mathcal{F} \mathcal{S}_{\bm{\mathrm{A}}^{-1}}\bed\\
&=
  \frac{1}{ | \det \bm{\mathrm{A}} | } \frac{1}{|  \det \bm{\mathrm{A}}^{-1} | } \mathcal{S}_{\bm{\mathrm{A}}^\Tr}  \mathcal{F} \bed \\
&= | \det \bm{\mathrm{A}}^{-1} | \frac{\mathcal{S}_{\bm{\mathrm{A}}^\Tr}}{ | \det \bm{\mathrm{A}}^{-1} |}   \bed \\
&=
  | \det \bm{\mathrm{A}}^{-1} | \, \bed_{\bm{\mathrm{A}}^{-\Tr} \mathbb{Z}^d}
\end{aligned}
$$

Note that the Fourier transform of a bed of impulses on the lattice
$\bm{\mathrm{A}}\mathbb{Z}^d$ is proportional to a bed of impulses on
the lattice $\bm{\mathrm{A}}^{-\Tr}\mathbb{Z}^d$. Note that since
$(\bm{\mathrm{A}}^{-\Tr})^\Tr \bm{\mathrm{A}} = \bm{\mathrm{A}}^{-1} \bm{\mathrm{A}} = \bm{\mathrm{I}}$,
the columns of ${\bm{\mathrm{A}}}$ and ${\bm{\mathrm{A}}}^{-\Tr}$ are
orthonormal. It is convenient to define this particular lattice as the
*dual* or *reciprocal* of $\mathcal{L}$,

$$
\mathcal{L}^\ast = \bm{\mathrm{A}}^{-\Tr} \mathbb{Z}^d.
$$ 

Note that

$$
\mathrm{vol}\, \mathcal{L}^\ast = |\det {\bm{\mathrm{A}}}^{-1} | = \frac{1}{\mathrm{vol}\, \mathcal{L}}.
$$

In this notation, we can write 
$$
\begin{aligned}
 \bed_{\mathcal{L}}  \fff  \mathrm{vol}\, \mathcal{L}^\ast \,  \bed_{\mathcal{L}^\ast}  
\end{aligned}
$$

### Sampling theorem on a lattice

Let ${\bm{\mathrm{A}}}$ be a regular matrix defining the lattice
$\mathcal{L} = \bm{\mathrm{A}} \mathbb{Z}^d$. We will denote by
$\mathcal{L}_0 = \bm{\mathrm{A}}\left[-\frac{1}{2},\frac{1}{2}\right]^d$
the *structural element* of the lattice, i.e., the smallest set such
that

$$
\bigcup_{\bm{\mathrm{p}} \in \mathcal{L} } (\bm{\mathrm{p}} + \mathcal{L}_0)  = \mathbb{R}^d.
$$

We will say that a function $f$ is *$\mathcal{L}_0$-band limited* if
$ \mathrm{supp}(F) \subset\mathcal{L}_0$. 

Let us denote
${\bm{\mathrm{B}}} = {\bm{\mathrm{A}}}^{-\Tr}$ and define
$g(\bm{\mathrm{x}}) = f(\bm{\mathrm{B}} \bm{\mathrm{x}}) = (\mathcal{S}_{\bm{\mathrm{B}}} f)(\bm{\mathrm{x}})$.

According to the stretching property of the Fourier transform,

$$
G(\bm{\mathrm{\xi}}) = \frac{1}{|\det \bm{\mathrm{B}}| } F(\bm{\mathrm{B}}^{-\Tr} \bm{\mathrm{\xi}})  = | \det \bm{\mathrm{A}} | \, F(\bm{\mathrm{A}}\bm{\mathrm{\xi}}).
$$

Hence,

$$
\mathrm{supp}(G) = \{ {\bm{\mathrm{\xi}}} : F({\bm{\mathrm{A}}}{\bm{\mathrm{\xi}}}) \ne 0  \} = {\bm{\mathrm{A}}}^{-1}  \{ {\bm{\mathrm{\xi}}} : F({\bm{\mathrm{\xi}}}) \ne 0  \} = {\bm{\mathrm{A}}}^{-1} \mathrm{supp}(F) \subset \left[-\frac{1}{2},\frac{1}{2}\right]^d.
$$

Hence, $g$ is band-limited on $\left[-\frac{1}{2},\frac{1}{2}\right]^d$
and we can write with
${\bm{\mathrm{y}}} = {\bm{\mathrm{B}}}{\bm{\mathrm{x}}}$

$$
\begin{aligned}
f(\bm{\mathrm{y}}) &= f(\bm{\mathrm{B}}\bm{\mathrm{x}}) = g(\bm{\mathrm{x}}) = \sum_{\bm{\mathrm{p}} \in \mathbb{Z}^d} g(\bm{\mathrm{p}}) \, \mathrm{sinc}(\bm{\mathrm{x}} - \bm{\mathrm{p}}) \\
&= \sum_{\bm{\mathrm{p}} \in \mathbb{Z}^d} g(\bm{\mathrm{p}}) \, \mathrm{sinc}(\bm{\mathrm{B}}^{-1} \bm{\mathrm{y}} - \bm{\mathrm{p}}) 
= \sum_{\bm{\mathrm{p}} \in \mathbb{Z}^d} g(\bm{\mathrm{p}}) \, \mathrm{sinc}(\bm{\mathrm{B}}^{-1} ( \bm{\mathrm{y}} - \bm{\mathrm{B}}\bm{\mathrm{p}})) \\
&= \sum_{\bm{\mathrm{q}} \in \bm{\mathrm{B}} \mathbb{Z}^d} g(\bm{\mathrm{B}}^{-1}\bm{\mathrm{q}}) \, \mathrm{sinc}(\bm{\mathrm{B}}^{-1} ( \bm{\mathrm{y}} - \bm{\mathrm{q}})) \\
&= \sum_{\bm{\mathrm{q}} \in \mathcal{L}^\ast } f(\bm{\mathrm{q}}) \, \mathrm{sinc}(\bm{\mathrm{A}}^{\Tr} ( \bm{\mathrm{y}} - \bm{\mathrm{q}}))
\end{aligned}
$$

In other words, a $\mathcal{L}_0$-band limited $f$ can be perfectly
reconstructed from its samples on $\mathcal{L}^\ast$.

### Anti-aliasing

Many real-world signals such as images are not truly band-limited
(though do exhibit decaying spectrum). In order to avoid aliasing at
sampling, the signal has to be made band-limited. If sampling occurs on
lattice $\mathcal{L}^\ast$, the signal has to be $\mathcal{L}$
band-limited in order to allow perfect reconstruction. A real signal is
made band-limited by applying a low pass *anti-aliasing filter*
immediately prior to sampling. The filter has to be applied in the
analog domain – in electronics when sampling a one-dimensional time
signal, or in optics when sampling an image. A lens has a natural high
frequency cut-off due to diffraction limit, which acts as an
anti-aliasing filter.

The ideal anti-aliasing filter should have the frequency response
$H_{\mathcal{L}}(\bm{\mathrm{\xi}})$ with $\mathrm{supp}(H_{\mathcal{L}}) \subset \mathcal{L_0}$, 
which is achieved by
$H_{\mathcal{L}}(\bm{\mathrm{\xi}}) = \mathrm{box}( \bm{\mathrm{A}}^{-1} \bm{\mathrm{\xi}} )$.
Using the stretching property of the Fourier transform, we obtain the
impulse response of the filter,

$$
h_{\mathcal{L}} (\bm{\mathrm{x}})  = \mathrm{vol}(\mathcal{L}) \, \mathrm{sinc}(\bm{\mathrm{A}}^\Tr \bm{\mathrm{x}}).
$$

Real anti-aliasing filters are never ideal, therefore some amount of
aliasing is always present in real systems.

## Interpolation


Sampling theorem tells us that when a function $f$ is band-limited on
the lattice $\mathcal{L} = {\bm{\mathrm{A}}} \mathbb{Z}^d$ (that is,
$ \mathrm{supp}(F) \subset {\bm{\mathrm{A}}}\left[-\frac{1}{2},\frac{1}{2}\right]^d$),
it can be perfectly reconstructed from its samples
$f|_{\mathcal{L}^\ast}$ on the dual lattice
$\mathcal{L}^\ast = {\bm{\mathrm{A}}}^{\ast} \mathbb{Z}^d$, where
${\bm{\mathrm{A}}}^\ast = {\bm{\mathrm{A}}}^{-\Tr}$, by the following
formula 

$$
\begin{aligned}
f(\bm{\mathrm{x}}) &=  \sum_{\bm{\mathrm{p}} \in \mathcal{L}^\ast } f(\bm{\mathrm{p}}) \, \mathrm{sinc}(\bm{\mathrm{A}}^{\Tr} ( \bm{\mathrm{x}} - \bm{\mathrm{p}})) 
=  \sum_{\bm{\mathrm{n}} \in \mathbb{Z}^d} f(\bm{\mathrm{A}}^\ast \bm{\mathrm{n}}) \, \mathrm{sinc}(\bm{\mathrm{A}}^{\Tr} \bm{\mathrm{x}} - \bm{\mathrm{n}}).
\end{aligned}
$$

Let us examine the expression

$$
\mathrm{sinc}({\bm{\mathrm{A}}}^{\Tr} {\bm{\mathrm{x}}} - {\bm{\mathrm{n}}}) = \mathrm{sinc}(\langle {\bm{\mathrm{a}}}_1, {\bm{\mathrm{x}}} \rangle - n_1) \cdots  \mathrm{sinc}(\langle {\bm{\mathrm{a}}}_d, {\bm{\mathrm{x}}} \rangle - n_d).
$$

Remembering that $\mathrm{sinc}(0) = 1$ and $\mathrm{sinc}(n) = 0$ for
every $n \in \mathbb{Z} \setminus \{ 0 \}$, we obtain that for every
${\bm{\mathrm{p}}} \in \mathcal{L}^\ast$, we can write
${\bm{\mathrm{p}}} = {\bm{\mathrm{A}}}^\ast {\bm{\mathrm{k}}}$ for
${\bm{\mathrm{k}}} \in \mathbb{Z}^d$, obtaining 

$$
\begin{aligned}
\mathrm{sinc}({\bm{\mathrm{A}}}^{\Tr} {\bm{\mathrm{A}}}^\ast {\bm{\mathrm{k}}} - {\bm{\mathrm{n}}}) &= \mathrm{sinc}({\bm{\mathrm{k}}} - {\bm{\mathrm{n}}}) = \mathrm{sinc}(k_1 - n_1) \cdots  \mathrm{sinc}(k_d - n_d) \\
&= \left\{ \begin{array}{cl} 1 & : {\bm{\mathrm{k}}} = {\bm{\mathrm{n}}} \\ 0 & : \mathrm{otherwise.} \end{array}  \right.
\end{aligned}
$$

This means that on the points of the lattice $\mathcal{L}^\ast$, the
interpolation formula yields exactly the sampled values.

### Interpolation kernels

The problem with the sinc is that it has infinite support. In real
systems, the interpolation formula takes place in an analog piece of
hardware (e.g., electronics or optics), which is limited to simpler
functions with finite support. A real interpolation formula is therefore
a revisited version of its ideal counterpart,

$$
f(\bm{\mathrm{x}}) =    \sum_{\bm{\mathrm{n}} \in \mathbb{Z}^d} f(\bm{\mathrm{A}}^\ast \bm{\mathrm{n}}) \, p(\bm{\mathrm{A}}^{\Tr} \bm{\mathrm{x}} - \bm{\mathrm{n}}),
$$

where $p$ is a $d$-dimensional *interpolation* or *reconstruction
kernel|* having the property that $p(\bm{\mathrm{0}}) = 1$ and
$p|_{\mathbb{Z}^d \setminus \{  \bm{\mathrm{0}}\} } = 0$. In practice,
a *separable* kernel $p(\bm{\mathrm{x}}) = p(x_1)\cdots p(x_d)$ is
often used. The one-dimensional kernel $p$ again has to satisfy
$p(0) = 1$ and $p(n)=0$ for $n \in \mathbb{Z} \setminus \{0\}$. The
following interpolation kernels are frequently used:

-   *Zeroth-order* or *nearest-neighbor interpolation*:
    $p(x) = \mathrm{rect}(x)$. In case of one-dimensional signals, a
    causal version of this kernel,
    $p(x) = \mathrm{rect}\left(x-\frac{1}{2}\right)$ is used.
    Zeroth-order interpolation sets the value of the continuous-domain
    function to the value of the nearest (causal in the latter case)
    discrete sample.

-   *First-order* or *(multi-) linear interpolation*:
    $p(x) =  \max (  1 - |x|, 0 )$. First-order interpolation continues
    the function between the samples linearly, essentially connecting
    the samples by lines (hyperplanes in the $d$-dimensional case).
    Bi-linear interpolation is the particular case for $d=2$.

-   *Third-order* or *(multi-) cubic interpolation*:

    $$
    p(x) = \left\{ \begin{array}{ll} 
    (a+2) |x|^3 - (a+3) |x|^2 + 1  & : |x| \le 1 \\
    a|x|^3 - 5a|x|^2 + 8a|x| - 4a & : 1 < |x| < 2 \\
    0 & : \mathrm{otherwise}.
    \end{array}  \right.
    $$ 

    The parameter $a$ is usually set to $a=-0.5$
    or $a=-0.75$. Cubic interpolation essentially fits a cubic Hermite
    spline to the data. Bi-cubic interpolation is the particular case
    for $d=2$.

-   *Lanczos interpolation*: 
	
	$$
	p(x) = \left\{ \begin{array}{ll} 
    \mathrm{sinc}(x) \, \mathrm{sinc}\left( \frac{x}{a} \right)  & : |x| < a \\
    0 & : \mathrm{otherwise}.
    \end{array}  \right.
    $$

    The parameter $a$ is a positive integer,
    typically 2 or 3, which determines the size of the kernel. The
    Lanczos kernel has $2a − 1$ lobes: a positive one at the center, and
    $a − 1$ alternating negative and positive lobes on each side. As
    long as the parameter $a$ is a positive integer, the Lanczos kernel
    is continuous everywhere, and its derivative is defined and
    continuous everywhere (including $x = \pm a$, where both sinc
    functions go to zero). Therefore, the reconstructed signal too will
    be $\mathcal{C}^1$.

[^1]: Technically, test functions must be sufficiently well-behaved. We
    will not dwell on these subtleties; in practice, any smooth
    compactly supported function will suffice for our purpose.

[^2]: In finite-dimensional spaces, the space of linear functionals is
    isomorphic to the vector space itself, and distributions are just a
    funny way to define vectors through their inner products with test
    vectors. On the other hand, in infinitely-dimensional spaces, such
    as our spaces of functions, no such isomorphism exists.
    Consequently, every function can be described as a distribution, but
    there are distributions such as the $\delta$ that cannot be
    described by a function. It can be approximated to an desired
    accuracy by a sequence of functions, but the limit point of the
    sequence is not a function.
