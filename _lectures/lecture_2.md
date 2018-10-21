---
title: "Lecture 1: Multidimensional signals and systems"
excerpt: Notation, Multidimensional signals and systems, Fourier transform
author: alex
---

## Notation 

In our course, we will deal almost exclusively with real functions. 
A scalar function will be denoted as $f : \RR^d \rightarrow \RR$, $f(\bb{x})$,
or simply $f$, with $d$ denoting the domain dimension. 
$d$-dimensional multi-indices will be denoted by $\bb{n} = (n_1,\dots,n_d) \in \bb{Z}^d$.
An operator acting on a function will be denoted by $\mathcal{H}(f)$ or simply $\mathcal{H}f$.


## Multidimensional signals and systems


A $d$-dimensional \emph{signal} is a real- (or sometimes complex-) valued function $f : \RR^d \rightarrow \RR$. In the case of $d=2$, 
we will refer to such signals as ''images'' and interpret the argument of the function, $\bb{x} = (x_1,x_2)$ as 
the spatial location of a point in an image; for $d=3$, we will refer to such signals as to images, interpreting 
the additional third dimension as ''time''. In the one-dimensional case,  the scalar argument will be denoted by $t$ and interpreted as ''time''. 

Denoting the space of real-valued signals as $\mathbb{F}(\RR^d,\RR) = \{f  :  f : \RR^d \rightarrow \RR \}$, a \emph{system} is an operator $\mathcal{H} : \mathbb{F}(\RR^d,\RR) \rightarrow \mathbb{F}(\RR^d,\RR)$. That is, a system $\mathcal{H}$ accepts a signal  as the input and produces another signal as the output.
%
For example the system receiving $f(\bb{x})$ and producing $f(\bb{x}-\bb{p})$ will be called \emph{translation} (or \emph{shift}) and the corresponding operator will be denoted as $\tau_{\bb{p}}$. With this notation we can write $(\tau_{\bb{p}} f)(\bb{x}) = f(\bb{x}-\bb{p})$. 

### Linearity

A system $\mathcal{H}$ is said to be a \emph{linear system} is it is additive and homogenous, i.e., for every $f,g \in \mathbb{F}(\RR^d,\RR)$ and $a, b \in \RR$,

$$
\mathcal{H}(a f + b g) = a \mathcal{H}f + b \mathcal{H}g .
$$

The output of (almost) every linear system can be represented as a linear functional of the form

$$
(\mathcal{H}f)(\bb{x}) = \int_{\RR^d} h(\bb{x},\bb{x}') f(\bb{x}') d\bb{x}',
$$

where the function (more generally, distribution -- see in the sequel)  $h(\bb{x},\bb{x}')$ is called the \emph{impulse response} of the system. Informally, $h(\bb{x},\bb{x}')$ tells what will be the output of the system at point $\bb{x}$ when the input is an impulse at $\bb{x}'$.  
%
An alternative way to view the above description is by defining the standard inner product on $L^2(\RR) \subset \mathbb{F}(\RR^d,\RR)$:

$$
\langle f,g \rangle = \int_{\RR^d} f(\bb{x}) g(\bb{x})^\ast d\bb{x};
$$

since most of the time we will be interested in real-valued functions, we will often omit the complex conjugate from the second argument.

### Shift invariance

A linear system $\mathcal{H}$ is called \emph{shift invariant} (LSI for short) if it commutes with the translation operator, i.e., for every $\bb{p} \in \RR^d$

$$
\tau_{\bb{p}} \mathcal{H} = \mathcal{H}\tau_{\bb{p}}.
$$

In other words, the output of the system given a shifted input is the same as the output of the system shifted by the same amount. 
Substituting a specific input $f$ and shift $\bb{p}$, we have the identity


<!-- \begin{eqnarray*} -->
$$
\int_{\RR^d} h(\bb{x} - \bb{p},\bb{x}') f(\bb{x}') d\bb{x}' = (\mathcal{H}f (\bb{x}) )(\bb{x}-\bb{p}) = (\mathcal{H}f (\bb{x} - \bb{p}) )(\bb{x}) = \int_{\RR^d} h(\bb{x},\bb{x}') f(\bb{x}'- \bb{p}) d\bb{x}' \\
= \int_{\RR^d} h(\bb{x},\bb{x}'' + \bb{p}) f(\bb{x}'') d\bb{x}''.
$$
<!-- \end{eqnarray*} -->


Since the latter holds for every $h$ and $\bb{p}$, we have $h(\bb{x}-\bb{p}, \bb{x}') = h(\bb{x}, \bb{x}'+\bb{p})$ at every $\bb{x}$. In other words, $h(\bb{x},\bb{x}')$ is effectively only a function of $\bb{x}-\bb{x}'$, which we will continue denoting as ''h'' with some abuse of notation. This particular structure is called Toeplitz and is the multidimensional infinite support equivalent of a circulant matrix.  The response of an LSI system can be therefore represented as


$$
(\mathcal{H}f)(\bb{x}) = \int_{\RR^d} h(\bb{x}-\bb{x}') f(\bb{x}') d\bb{x}'.
$$


The latter linear operation is called \emph{convolution} and denotes as $h \ast f$. It is straightforward to show that $h \ast f = f \ast h$. 
For real-valued functions, convolution of $h$ with $f$ can be also seen as $(h \ast f)(\bb{x}) = \langle \tau_{\bb{x}} \bar{h}, f  \rangle$,
where $\bar{h}(\bb{x}) = h(-\bb{x})$.

## Fourier transform

### Harmonics

A multi-dimensional \emph{harmonic} is the function $\phi_{\bb{\xi}}(\bb{x}) = e^{2\pi \ii \, \bb{\xi}^\Tr \bb{x}}$ with the $d$-dimensional real parameter $\bb{\xi}$. 
Note that the function has constant values along lines perpendicular to the direction of $\bb{\xi}$, and makes a full period every $1/\| \bb{\xi} \|$ units of space in the direction of $\bb{xi}$. In other words,  $\phi_{\bb{\xi}}$ looks like a periodic wave in the direction of $\bb{\xi}$ with the wavelength inversely proportional to $\|\bb{\xi}\|$ (or frequency proportional to $\|\bb{\xi}\|$). For this reason, $\bb{\xi}$ is usually referred to as the \emph{spatial frequency} measured in inverse length units.  
Note that $\phi_{\bb{\xi}}(\bb{p}) = \phi_{\bb{p}}(\bb{\xi})$, meaning that we can freely exchange the roles of the space location $\bb{x}$ and the spatial frequency $\bb{\xi}$. 


### Fourier transform

Orthogonal projections onto the set of harmonics is a function of the argument $\bb{\xi}$,

$$
F(\bb{\xi}) = \langle f, \phi_{\bb{\xi}} \rangle = \int_{\RR^d} f(\bb{x}) e^{-2\pi \ii \, \bb{\xi}^\Tr \bb{x}} d\bb{x},
$$

is called the (forward) \emph{Fourier transform} of $f$. We will denote it as the action of the operator $\mathcal{F}$ on $f$ or simply by the capital $F$. The argument of $F$  is conventionally referred to as \emph{frequency}, while the argument of $f$ is called \emph{space}. 

The inverse Fourier transform is defined as

$$
f(\bb{x}) =  (\mathcal{F}^{-1} F)(\bb{x}) =  \int_{\RR^d} F(\bb{\xi}) e^{2\pi \ii \, \bb{\xi}^\Tr \bb{x}} d\bb{\xi}.
$$

In what follows, we prove several very useful properties of the Fourier transforms.


### Tensor products

Let $f_1,\dots, f_d : \RR \rightarrow \RR$ be functions of a single variable, and let us define the tensor product $(f_1 \otimes \cdots \otimes f_d)(\bb{x}) = f_1(x_1) \cdots f_d(x_d)$. Applying the Fourier transform yields

\begin{eqnarray*}
(\mathcal{F}(f_1 \otimes \cdots \otimes f_d))( \bb{\xi} ) &=& \int_{\RR^d} f_1(x_1) \cdots f_d(x_d) e^{-2\pi \ii (x_1 \xi_1 + \cdots + x_d \xi_d) } dx_1 \cdots dx_d \nonumber\\
&=&  \int_\RR f_1(x_1) e^{-2\pi \ii x_1 \xi_1} dx_1 \cdots  \int_\RR f_d(x_d) e^{-2\pi \ii x_d \xi_d} dx_d \nonumber\\
&=&(\mathcal{F} f_1) (\xi_1) \cdots (\mathcal{F} f_d) (\xi_d).
\end{eqnarray*}

We can write this result in the following convenient notation:

\begin{align*}
\Aboxed{  f_1 \otimes \cdots \otimes f_d & \fff  F_1 \cdots F_d  }
\end{align*}

A useful example is the $d$-dimensional box function that can be defined as a tensor product of one-dimensional rectangular functions:
$\mathrm{box} = \mathrm{rect} \otimes \cdots \otimes \mathrm{rect}$, where

$$
\mathrm{rect}(x) = \left\{ \begin{array}{lcl} 0 & : & |x| > \frac{1}{2} \\
1 & : & |x| < \frac{1}{2}.
\end{array} \right.
$$

Since $(\mathcal{F} \mathrm{rect})(\xi) = \mathrm{sinc}(\xi) = \frac{\sin(\pi \xi)}{\pi \xi} $, we have $(\mathcal{F} \mathrm{box})(\bb{\xi}) =  \mathrm{sinc}(\xi_1) \cdots \mathrm{sinc}(\xi_d)$.


### Translation and modulation

 Let $\bb{p} \in \RR^d$ and $f : \RR^d \rightarrow \RR$. Then,

\begin{eqnarray*}
(\mathcal{F}( \tau_{\bb{p}} f))( \bb{\xi} ) &=& \int_{\RR^d}  f(\bb{x}-\bb{p}) e^{-2\pi \ii \, \bb{\xi}^\Tr \bb{x}} d\bb{x} =  \int_{\RR^d}  f(\bb{y}) e^{-2\pi \ii \, \bb{\xi}^\Tr (\bb{y}+\bb{p})} d\bb{y} \\
&=& e^{-2\pi \ii \, \bb{\xi}^\Tr \bb{p}} \int_{\RR^d}  f(\bb{y}) e^{-2\pi \ii \, \bb{\xi}^\Tr \bb{y}} d\bb{y} = e^{-2\pi \ii \, \bb{\xi}^\Tr \bb{p}} (\mathcal{F}f)( \bb{\xi} )
\end{eqnarray*}

We can write this result as the duality between translation and \emph{modulation} (i.e., multiplication by a harmonic)

\begin{align*}
\Aboxed{ \tau_{\bb{p}} f & \fff  \phi_{-\bb{p}} F   }
\end{align*}

Obviously, the same relation holds when swapping the space and frequency domains.
%$\phi_{\bb{p}} f \fff  \tau_{\bb{p}} F$.
In other words, shift in the space domain results in modulation (the addition of linear phase) in the frequency domain and modulation in the space domain results in shift in the frequency domain. 

### Convolution

Let $\bb{p} \in \RR^d$ and $f,h : \RR^d \rightarrow \RR$. Then,

\begin{eqnarray*}
(\mathcal{F}( h \ast f))( \bb{\xi} ) &=& \int_{\RR^d} \left(  \int_{\RR^d}  h(\bb{x}-\bb{y}) f(\bb{y}) d\bb{y} \right) e^{-2\pi \ii \, \bb{\xi}^\Tr \bb{x}} d\bb{x}  \\
&=& \int_{\RR^d}  \int_{\RR^d}  h(\bb{x}') f(\bb{y})  e^{-2\pi \ii \, \bb{\xi}^\Tr (\bb{x}' + \bb{y})} d\bb{x}' d\bb{y} \\
&=&  \int_{\RR^d} h(\bb{x}')   e^{-2\pi \ii \, \bb{\xi}^\Tr \bb{x}'} d\bb{x}'  \cdot \int_{\RR^d}  f(\bb{y}) e^{-2\pi \ii \, \bb{\xi}^\Tr  \bb{y}} d\bb{y}  = 
(\mathcal{F} h )( \bb{\xi} )  \cdot (\mathcal{F} f )( \bb{\xi} ) 
\end{eqnarray*}

In other words, convolution is dual to multiplication

\begin{align*}
\Aboxed{h \ast f & \fff  H \cdot F   }
\end{align*}

Obviously, the relation $h \cdot f \fff H \ast F$ also holds.

The fact that convolution (a Toeplitz operator) becomes a pointwise product (a diagonal operator) under the Fourier transform is the manifestation of the fact that the harmonics form an eigenbasis for every linear shift invariant system (or, equivalently, the Fourier transform diagonalizes every linear shift invariant operator). To see that, observe that when a harmonic $\phi_{\bb{\xi}}$ is given as the input to a system described by the impulse response $h$, the output,


\begin{eqnarray*}
(h \ast \phi_{\bb{\xi}}) (\bb{x}) &=& (\phi_{\bb{\xi}} \ast h) (\bb{x}) = \int_{\RR^d}   \phi_{\bb{\xi}}( \bb{x} - \bb{x}' )  h(\bb{x}')  d\bb{x}' = 
\int_{\RR^d}  h(\bb{x}') e^{2\pi \ii \, \bb{\xi}^\Tr  (\bb{x}-\bb{x}')}  d\bb{x}' \\
&=& \int_{\RR^d}  h(\bb{x}') e^{-2\pi \ii \, \bb{\xi}^\Tr \bb{x}'}  d\bb{x}' \cdot  e^{2\pi \ii \, \bb{\xi}^\Tr  \bb{x} } = F(\bb{\xi}) \phi_{\bb{\xi}} (\bb{x}),
\end{eqnarray*}


meaning that $\phi_{\bb{\xi}}$ is an eigenfunction of the system corresponding to the eigenvalue $F(\bb{\xi})$. Note that the spectrum is continuous and depends on the specific $h$, while the eigenfunctions are independent of $h$.\footnote{It is a general fact from linear algebra that operators can be jointly diagonalized if and only if they commute with each other.} 

### Stretching 

 Let $\bb{A} \in \RR^{d \times d}$ be a regular (i.e., non-singular) matrix and let us define the \emph{stretching} operator $\mathcal{S}_{\bb{A}} : f(\bb{x}) \mapsto f(\bb{A}\bb{x})$. Then,

\begin{eqnarray*}
(\mathcal{F}(  \mathcal{S}_{\bb{A}}  f))( \bb{\xi} ) &=& \int_{\RR^d}  f(\bb{A}\bb{x}) e^{-2\pi \ii \, \bb{\xi}^\Tr \bb{x}} d\bb{x} = 
\int_{\RR^d}  f(\bb{y}) e^{-2\pi \ii \, \bb{\xi}^\Tr \bb{A}^{-1} \bb{y}} \frac{d\bb{y} }{| \det \bb{A} |} = \frac{(\mathcal{F}f)(\bb{A}^{-\Tr} \bb{\xi} ) }{| \det \bb{A} |}.
\end{eqnarray*}


In other words, 

\begin{align*}
\Aboxed{\mathcal{S}_{\bb{A}} f & \fff \frac{\mathcal{S}_{\bb{A}^{-\Tr }} F }{ | \det \bb{A} | }     }
\end{align*}

This identity will reappear when talking about dual lattices in relation to sampling.

A useful example of the stretching property is the multi-dimensional Gaussian. Let us first remind ourselves (or accept without a proof) that in the one-dimensional case $e^{-x^2} \fff e^{-\pi \xi^2}$. The tensor product property immediately yields $e^{-\bb{x}^\Tr \bb{x} } \fff e^{-\pi \bb{\xi}^\Tr \bb{\xi} }$ in the case of unit covariance (a.k.a. normal) multidimensional Gaussian. A general Gaussian $e^{-\bb{x}^\Tr \bb{C}^{-1} \bb{x} }$ with the covariance matrix 
$\bb{C}$ can be obtained by stretching the normal Gaussian with the symmetric inverse square root matrix $\bb{A} = \bb{C}^{-\frac{1}{2}}$, since $(\bb{A} \bb{x})^\Tr (\bb{A} \bb{x}) = \bb{x}^\Tr \bb{A}^\Tr \bb{A} \bb{x} = \bb{x}^\Tr \bb{C}^{-1} \bb{x}$. The stretching property of the Fourier transform tells us that the corresponding stretch in frequency is by $\bb{A}^{-\Tr} = \bb{C}^{\frac{1}{2}}$ with furher scaling of the transform by the determinant of $\bb{C}^{\frac{1}{2}}$:

\begin{align*}
 e^{-\bb{x}^\Tr \bb{C}^{-1} \bb{x}} & \fff  \frac{1}{\sqrt{ \det \bb{C} } } e^{-\pi \bb{\xi}^\Tr \bb{C} \bb{\xi} }. 
\end{align*}

The covariance matrix $\bb{C}$ tells us how much the signal is concentrated in various spatial directions; its determinant can be used to quantify this spread: when $\det \bb{C}$ is small, the signal is localized in space. In the frequency domain, the signal remains Gaussian with the inverse covariance. Since both $\bb{C}$ and $\bb{C}^{-1}$ have the same eigenbasis but reciprocal eigenvalues, observe that the more the signal is concentrated in space in a certain direction, the more it is spread in frequency in the corresponding direction, and vice versa. The quantity $\pi  \det \bb{C}^{-1}$ measures the spread of the signal in frequency. Since $\det \bb{C} \cdot \pi  \det \bb{C}^{-1} = \pi$, we conclude that the Gaussian cannot be simultaneously localized both in space and frequency. This conclusion holds for every signal\footnote{Actually, the Gaussian achieves the best possible simultaneous localization in both domains.} -- a result known as the \emph{uncertainty principle} with profound implications on why the Universe looks like it looks.  


### Rotation

 An important particular case of the stretching property is when $\bb{A}$ is an orthonormal matrix representing rotations (and, possibly, reflections). Let us define the \emph{rotation} operator as $\mathcal{R}_{\bb{R}} : f(\bb{x}) \mapsto f(\bb{R}\bb{x})$, where $\bb{R}$ is a rotation matrix. Since orthonormal matrices satisfy $\bb{R}^{-1} = \bb{R}^\Tr$ and $ \det \bb{R} = \pm 1$, the stretching property reduces to $\mathcal{R}_{\bb{R}} f \fff \mathcal{R}_{\bb{R}} F$. In other words, the Fourier transform commutes with rotation:

\begin{align*}
\Aboxed{ \mathcal{F} \mathcal{R}_{\bb{R}}   = \mathcal{R}_{\bb{R}} \mathcal{F}  }.
\end{align*}
  
### Projection

 Let us define the \emph{projection} operator along the last spatial coordinate as 

$$
\mathcal{P} : f(\bb{x}) \mapsto \int_{\RR} f(x_1,\cdots,x_{d-1}, x_d) dx_d. 
$$

Note that projection maps functions from $\mathbb{F}(\RR^d,\RR)$ to 
$\mathbb{F}(\RR^{d-1},\RR)$. Interpreting $f$ as a $d$-dimensional probability density function, the latter operation can be interpreted as \emph{marginalization} with respect to $x_d$. 
%
Invoking the $(d-1)$-dimensional Fourier transform on $\mathcal{P} f$ yields

\begin{eqnarray*}
(\mathcal{F}( \mathcal{P} f))( \bb{\xi} ) &=& \int_{\RR^{d-1}}  \left(   \int_{\RR} f(x_1,\cdots,x_{d-1}, x_d) dx_d  \right)  e^{-2\pi \ii \, (x_1 \xi_1 + \cdots x_{d-1}\xi_{d-1} )} dx_1 \cdots dx_{d-1} \\
&=&\left. \int_{\RR^d} f(\bb{x})  e^{-2\pi \ii \, \bb{x}^\Tr \bb{\xi} } d\bb{x} \right|_{\xi_d = 0} = (\mathcal{F}f)(\xi_1,\dots, \xi_{d-1}, 0).
\end{eqnarray*}

Defining the \emph{slice} operator $Q : f(\bb{x}) \mapsto f(x_1,\dots, x_{d-1},0)$ from $\mathbb{F}(\RR^d,\RR)$ to $\mathbb{F}(\RR^{d-1},\RR)$ we can express the latter result more compactly as $\mathcal{F}\mathcal{P} = \mathcal{Q}\mathcal{F}$. Note that while on the right hand side $\mathcal{F}$ denotes a $d$-dimensional Fourier transform, on the left hand side it stands for the $(d-1)$-dimensional counterpart. 

Applying a rotation operator on the right to both sides of the identity yields
$\mathcal{F}\mathcal{P} \mathcal{R}_{\bb{R}} = \mathcal{Q}\mathcal{F} \mathcal{R}_{\bb{R}} = \mathcal{Q} \mathcal{R}_{\bb{R}} \mathcal{F}$, 

where in the last passage we used the commutativity of the Fourier transform with rotation. We interpret the composition  $\mathcal{P} \mathcal{R}_{\bb{R}}$ as a general projection operator,  $\mathcal{P}_{\bb{R}}$, that first rotates the function and then project along the last axis. This essentially allows to project the function along any direction. In the same manner, we interpret $\mathcal{Q} \mathcal{R}_{\bb{R}}$ as a general slice operator, $\mathcal{Q}_{\bb{R}}$, slicing the function along an arbitrary direction. 
This general result is known as the \emph{slice-projection theorem} that in our notation can be expressed as

\begin{align*}
\Aboxed{ \mathcal{F} \mathcal{P}_{\bb{R}}   = \mathcal{Q}_{\bb{R}} \mathcal{F}  }.
\end{align*}

An extremely important example where this result is used is computerized tomography (CT). In an idealized scenario, let us assume a $d=2$ dimensional world, in which we are interested in measuring the density of a slice of the human body, denoted by the function $f(x,y)$. Being unable to actually slice it (without going to jail), the next best thing we can do is to irradiate it with penetrating radiation (x-rays) from the side. Let us assume that an x-ray source sends parallel rays from one side of the body to the other side along the vertical ($y$) direction, where a linear detector measures the intensity profile $I(x)$. According to the Beer-Lambert law of light attenuation, the measured intensity is given by

$$
I(x) = I_0\, \exp\left(-\int_\RR f(x,y) dy \right),
$$

where $I_0$ is the emitted intensity. Taking the logarithm yields

$$
p(x) = -\log \frac{I(x)}{I_0} = \int_\RR f(x,y) dy = \mathcal{P} f.
$$

Rotating the emitter-detector setup around the body yields a collection of projections $\mathcal{P}_\theta f$ (note that in two dimensions, the rotation matrix is parameterized by a single angle $\theta$). The function $(x,\theta) \mapsto (\mathcal{P}_\theta f)(x)$ is often referred to as the \emph{Radon transform} or the \emph{sinogram} of $f$. The slice projection theorem tells us that the one-dimensional Fourier transform of each such projection $\mathcal{P}_\theta f$ yields a correspondingly directed slice $\mathcal{Q}_\theta f$ of the two-dimension Fourier transform of the unknown function $f$. Collecting enough projections, it is ''just'' a matter of numerics to estimate the said transform and invert it, yielding what we see on the screen as a slice of a CT scan. 
