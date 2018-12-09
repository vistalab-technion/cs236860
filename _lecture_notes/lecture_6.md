---
title: "Lecture 6: Patch-based priors"
excerpt: Patch-based denoising, NLM, Bilateral filter
author: alex
copyright: alex
---


We continue our quest for a better model for the prior probability
$P_ {\mathpzc{F}}$ to be used in MAP and Bayesian estimators. While it is
almost hopeless to be able to formulate such a prior on the entire
image, it is a much more practical task to do it on a small region of an
image (known in the image processing jargon as a *patch*). Let us first
fix some size $T$ and define a square domain
$\square = \left[-\frac{T}{2},\frac{T}{2}\right]^d$. A patch around some
point ${\bm{\mathrm{x}}}$ in an image $f$ can be thought of as
$\tau_ {\bm{\mathrm{x}}} f |_ \square$. We can therefore define a *patch
operator*
$\pi_ {\bm{\mathrm{p}}} : \mathbb{F}( {\bm{\mathrm{R}}}^d, {\bm{\mathrm{R}}}) \rightarrow \mathbb{F}( \square, {\bm{\mathrm{R}}})$
taking $f$ as the input and producing
$\tau_ {\bm{\mathrm{x}}} f |_ \square$ as the output. In these terms, we
redefine the problem of modeling the prior probability of a natural
image $\mathpzc{F}$ as the problem of modeling the prior probability of
$\pi_ {\bm{\mathrm{x}}} \mathpzc{F}$ for every location
${\bm{\mathrm{x}}}$ in the image domain. To that end, let us assume to
have access to a (potentially, very large) collection of clean patches
$\{p_ 1,\dots,p_ K\}$; such a collection can be obtained by taking a very
large representative set of images, decomposing them into a collection
of overlapping patches and clustering the latter into $K$ samples. We
will further assume that a patch in a natural image
$\pi_ {\bm{\mathrm{q}}} \mathpzc{F}$ is uniformly distributed over that
collection, i.e., $\pi_ {\bm{\mathrm{x}}} \mathpzc{F} = p_ k$ with
probability $\frac{1}{K}$ for $k=1,\dots,K$. Note that the simple
assumption of uniform distribution over the collection of example
patches leads to a potentially very intricate prior distribution of
$\pi_ {\bm{\mathrm{x}}} \mathpzc{F}$ itself, since the collection
captures the intricate structure of representative image patches!

## Patch-based denoising


### Denoising using multiple images

Let us simplify our archetypal inverse problem assuming that the latent
signal $\mathpzc{F}$ is only degraded by white additive noise,
$\mathpzc{Y} = \mathpzc{F} + \mathpzc{N}$, admitting a Gaussian
distribution
$\mathpzc{N}({\bm{\mathrm{x}}}) \sim \mathcal{N}(0,\sigma_ \mathpzc{N}^2)$.
This problem is known as *Gaussian denoising*. The likelihood of
$\pi_ {\bm{\mathrm{x}}} \mathpzc{Y}$ given that
$\pi_ {\bm{\mathrm{x}}} \mathpzc{F}$ was formed using the patch $p_ k$
is simply given by the density of $\pi_ \bm{\mathrm{x}} \mathpzc{N}$,

$$P_ {\pi_ {\bm{\mathrm{x}}} \mathpzc{Y}  | \pi_ {\bm{\mathrm{x}}} \mathpzc{F} } (\pi_ {\bm{\mathrm{x}}} \mathpzc{Y} = y | \pi_ {\bm{\mathrm{x}}} \mathpzc{F} = p_ k) = f_ {\mathpzc{N}}( y - p_ k )  \propto e^{-\frac{ \| y - p_ k  \|^2 }{2 \sigma_ \mathpzc{N}^2} },$$

where the norm is the $L^2$ norm on $\square$. Since we asserted uniform
distribution of $k$ over $\{1,\dots,K\}$, the Bayes formula readliy
yields


$$P_ {\pi_ {\bm{\mathrm{x}}} \mathpzc{F} |\pi_ {\bm{\mathrm{x}}} \mathpzc{Y}   } ( \pi_ {\bm{\mathrm{x}}} \mathpzc{F} = p_ k | \pi_ {\bm{\mathrm{x}}} \mathpzc{Y} = y ) = \frac{  e^{-\frac{ \| y - p_ k  \|^2 }{2 \sigma_ \mathpzc{N}^2} }   }{  \displaystyle{\sum_ {i=1}^K   e^{-\frac{ \| y - p_ i  \|^2 }{2 \sigma_ \mathpzc{N}^2} }   }} .$$


Using this posterior, we can define the MSE of the estimator $\hat{f}$
at point ${\bm{\mathrm{x}}}$ as

$$\epsilon( \hat{f}({\bm{\mathrm{x}}}) ) =  \mathbb{E} ( \hat{f}({\bm{\mathrm{x}}}) - p_ k({\bm{\mathrm{0}}}) )^2  =  \sum_ { k=1 }^K ( \hat{f}({\bm{\mathrm{x}}}) - p_ k({\bm{\mathrm{0}}}) )^2  P_ {\pi_ {\bm{\mathrm{x}}} \mathpzc{F} |\pi_ {\bm{\mathrm{x}}} \mathpzc{Y}   } ( \pi_ {\bm{\mathrm{x}}} \mathpzc{F} = p_ k | \pi_ {\bm{\mathrm{x}}} \mathpzc{Y}  =  \pi_ {\bm{\mathrm{x}}} y );$$


the MMSE estimator is therefore defined as


$$\hat{f}({\bm{\mathrm{x}}})   = \mathrm{arg}\min_ {\hat{f}({\bm{\mathrm{x}}})  } \epsilon( \hat{f}({\bm{\mathrm{x}}})  ) = \mathrm{arg} \min_ { \hat{f}({\bm{\mathrm{x}}})  }  \sum_ { k=1 }^K ( \hat{f}({\bm{\mathrm{x}}}) - p_ k({\bm{\mathrm{0}}}) )^2   h_ {k} ({\bm{\mathrm{x}}}  ),$$


where
$h_ {k}( {\bm{\mathrm{x}}}  ) =  P_ {\pi_ {\bm{\mathrm{x}}} \mathpzc{F} |\pi_ {\bm{\mathrm{x}}} \mathpzc{Y}   } ( \pi_ {\bm{\mathrm{x}}} \mathpzc{F} = p_ k | \pi_ {\bm{\mathrm{x}}} \mathpzc{Y}  =  \pi_ {\bm{\mathrm{x}}} y )$.
Note that the sum of the latter weights over $k$ is always one. Also
note that $ p_ k({\bm{\mathrm{0}}})$ gives the value at the center of the
patch. The latter problem has a simple closed-form solution as the
weighted average


$$\hat{f}({\bm{\mathrm{x}}}) = \sum_ { k=1 }^K  h_ {k} ( {\bm{\mathrm{x}}}   ) p_ k({\bm{\mathrm{0}}}) = \frac{   \displaystyle{\sum_ {k=1}^K e^{-\frac{ \| \pi_ {\bm{\mathrm{x}}} y - p_ k  \|^2 }{2 \sigma_ \mathpzc{N}^2}  }  p_ k({\bm{\mathrm{0}}}) }   }{  \displaystyle{\sum_ {k=1}^K   e^{-\frac{ \| \pi_ {\bm{\mathrm{x}}} y - p_ k  \|^2 }{2 \sigma_ \mathpzc{N}^2} }   }}$$


of the central values from exemplar patches $\{p_ 1,\dots,p_ K\}$, with
weights inversely proportional to the distance of the said patches from
the corresponding input patch centered at ${\bm{\mathrm{x}}}$. Since $K$
might be very big and the exponential anyway decays very fast in
$ \| \pi_ {\bm{\mathrm{x}}} y - p_ k  \|$, in practice the above sum is
approximated by taking a fixed number of (approximate) *nearest
neighbors)* in the sense of the inter-patch distance
$ \| \pi_ {\bm{\mathrm{x}}} y - p_ k  \|$.

### Non-local means

It is sometimes inconvenient to assume the availability of exemplar
patches $\{p_ k\}$ coming from an external collection of images; instead,
it is often desirable to use the patches from the image itself for the
purpose of defining the prior on $\pi_ {\bm{\mathrm{x}}} \mathpzc{F}$.
The former MMSE estimator can be readily formulated in this setup as


$$\hat{f}({\bm{\mathrm{x}}}) =  \frac{   \displaystyle{\int_ {\RR^d} e^{-\frac{ \| \pi_ {\bm{\mathrm{x}}} y -  \pi_ {\bm{\mathrm{x}}'} y \|^2 }{2 \sigma_ \mathpzc{N}^2}  } y({\bm{\mathrm{x}}}')  d{\bm{\mathrm{x}}}' }   }{  \displaystyle{\int_ {\RR^d}   e^{-\frac{ \| \pi_ {\bm{\mathrm{x}}} y - \pi_ {\bm{\mathrm{x}}'} y  \|^2 }{2 \sigma_ \mathpzc{N}^2} }  d{\bm{\mathrm{x}}}'  }}$$


or, more practically,


$$\hat{f}({\bm{\mathrm{x}}}) = \frac{   \displaystyle{\sum_ {k=1}^K e^{-\frac{ \| \pi_ {\bm{\mathrm{x}}} y - \pi_ {\bm{\mathrm{x}}_ k} y  \|^2 }{2 \sigma_ \mathpzc{N}^2}  }  y({\bm{\mathrm{x}}}_ k) }   }{  \displaystyle{\sum_ {k=1}^K   e^{-\frac{ \| \pi_ {\bm{\mathrm{x}}} y - \pi_ {\bm{\mathrm{x}}_ k} y  \|^2 }{2 \sigma_ \mathpzc{N}^2} }   }},$$


where
$\{ \pi_ {\bm{\mathrm{x}}_ 1} y, \dots, \pi_ {\bm{\mathrm{x}}_ K} y \}$
are the $K$ nearest neighbors of $\pi_ {\bm{\mathrm{x}}} y$ in the
noisy image itself $y$ (variants of the PatchMatch are very useful to
compute the nearest neighbors). Such a filter is known as *non-local
means* (NLM) – it works like a regular local averaging that supresses
noise (by the law of large numbers), yet, in the case of NLM averaging
is non-local. Oftentimes, an additional (very slowly decaying) spatial
weight is added to down-weigh spatially distant pixels,


$$\hat{f}({\bm{\mathrm{x}}}) = \frac{   \displaystyle{\sum_ {k=1}^K e^{-\frac{ \| \pi_ {\bm{\mathrm{x}}} y - \pi_ {\bm{\mathrm{x}}_ k} y  \|^2 }{2 \sigma_ \mathpzc{N}^2} - \frac{\| {\bm{\mathrm{x}}}-{\bm{\mathrm{x}}}' \|^2_ 2 }{\sigma^2_ \mathrm{s}}  }  y({\bm{\mathrm{x}}}_ k) }   }{  \displaystyle{\sum_ {k=1}^K   e^{-\frac{ \| \pi_ {\bm{\mathrm{x}}} y - \pi_ {\bm{\mathrm{x}}_ k} y  \|^2 }{2 \sigma_ \mathpzc{N}^2}  - \frac{\| {\bm{\mathrm{x}}}-{\bm{\mathrm{x}}}' \|^2_ 2 }{\sigma^2_ \mathrm{s}}   }   }},$$


where $\sigma^2_ \mathrm{s}$ controls the width of the spatial weight.

### Bilateral filter

A particularly renowned setting of NLM (that was actually devised prior
to NLM) is the case of point-wise patches, i.e.,
$\pi_ {\bm{\mathrm{x}}} : f \rightarrow f({\bm{\mathrm{x}}})$. In this
case, the NLM estimator simplifies to


$$\hat{f}({\bm{\mathrm{x}}}) =  \frac{   \displaystyle{\int_ {\RR^d} e^{-\frac{ (y({\bm{\mathrm{x}}}) - y({\bm{\mathrm{x}}}'))^2 }{2 \sigma_ \mathpzc{N}^2} - \frac{\| {\bm{\mathrm{x}}}-{\bm{\mathrm{x}}}' \|^2_ 2 }{\sigma^2_ \mathrm{s}}  }  y({\bm{\mathrm{x}}}')  d{\bm{\mathrm{x}}}' }   }{  \displaystyle{\int_ {\RR^d}   e^{-\frac{ ( y({\bm{\mathrm{x}}}) - y({\bm{\mathrm{x}}}') )^2 }{2 \sigma_ \mathpzc{N}^2} - \frac{\| {\bm{\mathrm{x}}}-{\bm{\mathrm{x}}}' \|^2_ 2 }{\sigma^2_ \mathrm{s}}  }  d{\bm{\mathrm{x}}}'  }},$$


which is known as the *bilateral filter*. Bilateral filter can be
thought of as a non-shift invariant linear filter,


$$\hat{f}({\bm{\mathrm{x}}}) =  \int_ {\RR^d} h({\bm{\mathrm{x}}},{\bm{\mathrm{x}}}')  y({\bm{\mathrm{x}}}')  d{\bm{\mathrm{x}}}'$$


with the spatially-varying unit-DC impulse response


$$h({\bm{\mathrm{x}}},{\bm{\mathrm{x}}}')  =  \frac{    e^{-\frac{ (y({\bm{\mathrm{x}}}) - y({\bm{\mathrm{x}}}'))^2 }{2 \sigma_ \mathpzc{N}^2} - \frac{\| {\bm{\mathrm{x}}}-{\bm{\mathrm{x}}}' \|^2_ 2 }{\sigma^2_ \mathrm{s}}  }   }{  \displaystyle{\int_ {\RR^d}   e^{-\frac{ ( y({\bm{\mathrm{x}}}) - y({\bm{\mathrm{x}}}') )^2 }{2 \sigma_ \mathpzc{N}^2} - \frac{\| {\bm{\mathrm{x}}}-{\bm{\mathrm{x}}}' \|^2_ 2 }{\sigma^2_ \mathrm{s}}  }  d{\bm{\mathrm{x}}}'  }}.$$


Another useful model is to consider the so-called *lifted* space
${\bm{\mathrm{R}}}^d \times {\bm{\mathrm{R}}}$, in which the graph of
the image $f$ resides. The graph can be modeled as a delta distribution
on ${\bm{\mathrm{R}}}^d \times {\bm{\mathrm{R}}}$ supported on
$({\bm{\mathrm{x}}},f({\bm{\mathrm{x}}})$ (note that every function $f$
has such a representation, but not every function or distribution on
${\bm{\mathrm{R}}}^d \times {\bm{\mathrm{R}}}$ represents a valid
function). Then, convolving the latter distribution with an LSI Gaussian
kernel,


$$h(({\bm{\mathrm{x}}},f)) \propto  \exp\left(- \left\|  \left( \frac{\bm{\mathrm{x}}}{ \sqrt{2} \sigma_ \mathpzc{N} }, \frac{f}{ \sqrt{2}  \sigma_ \mathrm{s} } \right) \right\|^2_ 2 \right),$$


and marginalizing over $f$ yields precisely the bilateral filter.

### NLM as MMSE

Note that in its original formulation, the NLM estimator is not an MMSE
estimator. This stems from the fact that when using the patches of the
image itself as exemplars to define the prior on
$\pi_ {\bm{\mathrm{x}}} \mathpzc{F}$, the patches are contaminated with
noise! Fixing a set of $K$ locations ${\bm{\mathrm{x}}}_ k$ from which
exemplar patches are taken, we now do not have access to $p_ k$ but
rather to their version contaminated by noise,
$q_ k = p_ k + \pi_ {\bm{\mathrm{x}}_ k} \mathpzc{N} = p_ k + \mathpzc{N}_ k$.

Let us repeat the derivation of the posterior probability distribution.
The likelihood of $\pi_ {\bm{\mathrm{x}}} \mathpzc{Y}$ given that
$\pi_ {\bm{\mathrm{x}}} \mathpzc{F}$ was formed using the patch $p_ k$
is simply given by the density of $\pi_ {\bm{\mathrm{x}}} \mathpzc{N}$,

$$P_ {\pi_ {\bm{\mathrm{x}}} \mathpzc{Y}  | \pi_ {\bm{\mathrm{x}}} \mathpzc{F} } (\pi_ {\bm{\mathrm{x}}} \mathpzc{Y} = y | \pi_ {\bm{\mathrm{x}}} \mathpzc{F} = p_ k) \propto e^{-\frac{ \| y - p_ k  \|^2 }{2 \sigma_ \mathpzc{N}^2} };$$

however, we no more have access to $p_ k$. Instead, we can write

$$\begin{aligned}
 \| y - p_ k  \|^2 &=&  \| y - q_ k + \mathpzc{N}_ k  \|^2 \approx \mathbb{E} \| y - q_ k + \mathpzc{N}_ k  \|^2 = \| y - q_ k\|^2 + \mathbb{E} \| \mathpzc{N}_ k \|^2 + 2 \mathbb{E} \langle y-q_ k,  \mathpzc{N}_ k \rangle\end{aligned}$$

By definition of random Gaussian noise, since the norm is computed on
$\square$ of volume $T^d$,
$\mathbb{E} \| \mathpzc{N}_ k \|^2 = T^d \sigma^2_ {\mathpzc{N}}$. Letting
$y = p_ k + \pi_ {\bm{\mathrm{x}}} \mathpzc{N}$, the third term can be
written as

$$\mathbb{E} \langle y-q_ k,  \mathpzc{N}_ k \rangle = \mathbb{E} \langle p_ k + \pi_ {\bm{\mathrm{x}
}} \mathpzc{N} - p_ k - \mathpzc{N}_ k,  \mathpzc{N}_ k \rangle =  \mathbb{E} \langle \pi_ {\bm{\mathrm{x}}} \mathpzc{N},  \mathpzc{N}_ k \rangle - \mathbb{E} \| \mathpzc{N}_ k \|^2.$$


Even though the denoised and the reference patches may overlap, we will
assume that $ \pi_ {\bm{\mathrm{x}}} \mathpzc{N}$ and $\mathpzc{N}_ k$
are statistically independent and, since they are zero mean, they are
orthogonal. This yields
$\mathbb{E} \langle y-q_ k,  \mathpzc{N}_ k \rangle = -T^d \sigma^2_ {\mathpzc{N}}$;
combining the previous results, we have

$$P_ {\pi_ {\bm{\mathrm{x}}} \mathpzc{Y}  | \pi_ {\bm{\mathrm{x}}} \mathpzc{F} } (\pi_ {\bm{\mathrm{x}}} \mathpzc{Y} = y | \pi_ {\bm{\mathrm{x}}} \mathpzc{F} = p_ k)  \propto e^{\frac{ - \| y - q_ k  \|^2 + T^d \sigma^2_ {\mathpzc{N}} }{2 \sigma_ \mathpzc{N}^2} }$$


Similarly to the noiseless exemplars case, the Bayes formula yields


$$P_ {\pi_ {\bm{\mathrm{x}}} \mathpzc{F} |\pi_ {\bm{\mathrm{x}}} \mathpzc{Y}   } ( \pi_ {\bm{\mathrm{x}}} \mathpzc{F} = p_ k | \pi_ {\bm{\mathrm{x}}} \mathpzc{Y} = y ) = \frac{  e^{-\frac{ \| y - q_ k  \|^2 + T^d \sigma^2_ {\mathpzc{N}}  }{2 \sigma_ \mathpzc{N}^2} }   }{  \displaystyle{\sum_ {i=1}^K   e^{-\frac{ \| y - q_ i  \|^2 + T^d \sigma^2_ {\mathpzc{N}}  }{2 \sigma_ \mathpzc{N}^2} }   }} = h_ k({\bm{\mathrm{x}}}).$$


Using this posterior, we again define the MSE of the estimator $\hat{f}$
at point ${\bm{\mathrm{x}}}$ as


$$\epsilon( \hat{f}({\bm{\mathrm{x}}}) ) =  \mathbb{E} ( \hat{f}({\bm{\mathrm{x}}}) - p_ k({\bm{\mathrm{0}}}) )^2  =  \sum_ { k=1 }^K ( \hat{f}({\bm{\mathrm{x}}}) - p_ k({\bm{\mathrm{0}}}) )^2  h_ k({\bm{\mathrm{x}}}).$$


In the previous case, this yielded an estimator in the form of the
weighted sum


$$\pi_ {\bm{\mathrm{x}}} \hat{f}  = \sum_ {k=1}^K c_ k  p_ k,$$ 

with $c_ k = h_ k({\bm{\mathrm{x}}})$ summing to $1$. However, note that $p_ k$
are now unaccessible, so we have to write


$$\pi_ {\bm{\mathrm{x}}} \hat{f}= \sum_ {k=1}^K c_ k  q_ k = \sum_ {k=1}^K c_ k  p_ k  + \sum_ {k=1}^K c_ k  \mathpzc{N}_ k,$$


with some other set of weights ${\bm{\mathrm{c}}}$ also summing to $1$.
This yields 

$$\begin{aligned}
\epsilon( {\bm{\mathrm{c}}} ) = & \sum_ { j=1 }^K \| \pi_ {\bm{\mathrm{x}}}\hat{f} - p_ k \|^2 h_ k({\bm{\mathrm{x}}})  = \sum_ { k=1 }^K \left\| \sum_ {i=1}^K c_ i  p_ i + \sum_ {i=1}^K c_ i  \mathpzc{N}_ i - p_ k\right\|^2  h_ k({\bm{\mathrm{x}}}) \\
\approx &  \sum_ { k=1 }^K {\bm{\mathrm{E}}} \left\| \sum_ {i=1}^K c_ i  p_ i + \sum_ {i=1}^K c_ i  \mathpzc{N}_ i - p_ k\right\|^2  h_ k({\bm{\mathrm{x}}}) \\
\approx &  \sum_ { k=1 }^K \left(  \left\| \sum_ {i=1}^K c_ i  p_ i  - p_ k\right\|^2  +  \sum_ {i=1}^K c^2_ i {\bm{\mathrm{E}}}\| \mathpzc{N}_ i \|^2
 \right) h_ k({\bm{\mathrm{x}}}) \\
 =& \sum_ { k=1 }^K  \left\| \sum_ {i=1}^K c_ i  p_ i  - p_ k\right\|^2   h_ k({\bm{\mathrm{x}}}) + T^d  \sigma_ \mathpzc{N}^2 {\bm{\mathrm{c}}}^\Tr{\bm{\mathrm{c}}} \\
 =& {\bm{\mathrm{c}}}^\Tr ( {\bm{\mathrm{P}}} + T^d  \sigma_ \mathpzc{N}^2 {\bm{\mathrm{I}}} ) {\bm{\mathrm{c}}} - 2{\bm{\mathrm{c}}}^\Tr {\bm{\mathrm{P}}} {\bm{\mathrm{h}}} + {\bm{\mathrm{w}}}^\Tr {\bm{\mathrm{h}}}, \end{aligned}$$


where ${\bm{\mathrm{P}}}$ is a $K \times K$ matrix with the elements
$({\bm{\mathrm{P}}})_ {ij} = \langle p_ i, p_ j \rangle$,
${\bm{\mathrm{h}}}$ is the $K$-dimensional vector with the elements
$h_ k({\bm{\mathrm{x}}})$, and ${\bm{\mathrm{w}}}$ is the $K$-dimensional
vector with the elements $({\bm{\mathrm{w}}})_ i = \|  p_ i \|^2$. Our
MMSE estimator is therefore given as the solution to the constrained
minimization problem


$${\bm{\mathrm{c}}}_ \ast = \mathrm{arg}\min_ {\bm{\mathrm{c}} } \epsilon( {\bm{\mathrm{c}}} )  \,\,\, 
\mathrm{s.t.}\,\,\, {\bm{\mathrm{c}}}^\Tr {\bm{\mathrm{1}}} = 1.$$


Defining the Lagrangian
$ \epsilon( {\bm{\mathrm{c}}} ) + \lambda {\bm{\mathrm{1}}}^\Tr {\bm{\mathrm{c}}}$,
differentiating w.r.t. ${\bm{\mathrm{c}}}$ and setting the gradient to
zero yields


$${\bm{\mathrm{c}}}_ \ast = ({\bm{\mathrm{P}}} + T^d  \sigma_ \mathpzc{N}^2 {\bm{\mathrm{I}}} )^{-1} ({\bm{\mathrm{P}}} {\bm{\mathrm{h}}} + \lambda {\bm{\mathrm{1}}})$$


with the constant


$$\lambda = \frac{1 -  {\bm{\mathrm{1}}}^\Tr ({\bm{\mathrm{P}}} + T^d  \sigma_ \mathpzc{N}^2 {\bm{\mathrm{I}}} )^{-1}  {\bm{\mathrm{P}}} {\bm{\mathrm{h}}} }{ {\bm{\mathrm{1}}}^\Tr ({\bm{\mathrm{P}}} + T^d  \sigma_ \mathpzc{N}^2 {\bm{\mathrm{I}}} )^{-1}  {\bm{\mathrm{1}}}}$$


satisfying the unit sum constraint.

Note the inversion of the potentially very big matrix
${\bm{\mathrm{P}}} + T^d  \sigma_ \mathpzc{N}^2 {\bm{\mathrm{I}}}$ makes
this solution impractical. However, assuming that for every $i$,
$\langle p_ i, p_ i \rangle = T^d s^2_ {\mathpzc{F}}$ and for every
$j \ne i$, $\langle p_ i, p_ j \rangle = \rho T^d s^2_ {\mathpzc{F}}$,
where
$s^2_ {\mathpzc{F}} = {\bm{\mathrm{E}}} \mathpzc{F}({\bm{\mathrm{x}}})$,
the matrix the matrix ${\bm{\mathrm{P}}}$ assumes the simple form
${\bm{\mathrm{P}}} = \alpha {\bm{\mathrm{I}}} + \beta {\bm{\mathrm{1}}}{\bm{\mathrm{1}}}^\Tr$,
where $\alpha =   T^d (1-\rho)s^2_ {\mathpzc{F}}$ and
$\beta = T^d \rho  s^2_ {\mathpzc{F}}$. The matrix to invert therefore
becomes


$${\bm{\mathrm{P}}} + T^d  \sigma_ \mathpzc{N}^2 {\bm{\mathrm{I}}}  = \gamma {\bm{\mathrm{I}}} + \beta {\bm{\mathrm{1}}}{\bm{\mathrm{1}}}^\Tr,$$


where $\gamma =  \alpha +  T^d \sigma_ \mathpzc{N}^2$. Using the
Sherman-Morrison matrix identity, the inverse can be expressed as


$$\begin{aligned}
( \gamma {\bm{\mathrm{I}}} + \beta {\bm{\mathrm{1}}}{\bm{\mathrm{1}}}^\Tr )^{-1} &=& \frac{1}{\gamma} \left( {\bm{\mathrm{I}}} - \frac{\beta}{\gamma + \beta K} {\bm{\mathrm{1}}}{\bm{\mathrm{1}}}^\Tr \right).\end{aligned}$$


Substituting the latter result to the formula for
${\bm{\mathrm{c}}}_ \ast$, after tedious and boring algebra, we arrive at


$${\bm{\mathrm{c}}}_ \ast = \frac{(1-\rho)s^2_ {\mathpzc{F}} }{ (1-\rho) s^2_ {\mathpzc{F}} + \sigma^2_ {\mathpzc{N}}} \, {\bm{\mathrm{h}}} 
+  \frac{\sigma^2_ {\mathpzc{N}} }{ (1-\rho) s^2_ {\mathpzc{F}} + \sigma^2_ {\mathpzc{N}}} \, \frac{1}{K} \, {\bm{\mathrm{1}}}.$$


Note that the weight vector is given by a weighted average (a convex
combination) of the posterior probabilities ${\bm{\mathrm{h}}}$ and the
uniform vector $\frac{1}{K} \, {\bm{\mathrm{1}}}$. Interpreting the
signal to noise ratio as


$$\mathrm{SNR} =  \frac{(1-\rho)s^2_ {\mathpzc{F}} }{ \sigma^2_ {\mathpzc{N}}},$$


we can rewrite


$$c_ k = \frac{ \mathrm{SNR} \cdot h_ k + \frac{1}{K} }{  \mathrm{SNR}  + 1}.$$


For high SNR ($\sigma^2_ {\mathpzc{N}}$ approaching $0$), we obtain the
classical NLM with $c_ k = h_ k$.
