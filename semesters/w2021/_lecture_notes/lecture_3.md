---
title: "Lecture 3: Discrete-domain signals and systems"
excerpt: Discrete-domain Fourier transform, Sampling, Aliasing, Interpolation, Polyphase representation 
author: alex
copyright: alex
---


## Discrete-domain Fourier transform

When a function $f$ is sampled on the lattice $\mathcal{L}^\ast$, its
Fourier transform is periodized on the reciprocal lattice $\mathcal{L}$.
Therefore, it is sufficient to describe the value of the transform on
the *structural element* of the lattice, which can be parametrized as
$\bb{\xi} = \bb{A} \bb{\omega}$ with
$\bb{\omega} \in \left[ -\frac{1}{2}, \frac{1}{2} \right]$. The vector
$\bb{\omega}$ is called the *normalized* or *digital frequency*.

Substituting the normalized frequency into the periodized and scaled
Fourier transform of $f$ yields 

$$\begin{aligned}
F[\bb{\omega}] =& \mathrm{vol}\, \mathcal{L}\cdot (\bed _{\mathcal{L}} \ast F)( \bb{A} \bb{\omega} ) =  \mathcal{F}\left( \bed _{\mathcal{L}^\ast} \cdot f   \right)( \bb{A} \bb{\omega} ) =   \sum _{ \bb{p} \in \mathcal{L}^\ast } f(\bb{p}) \mathcal{F} (\tau _{\bb{p}} \delta)( \bb{A} \bb{\omega} ) \\
=&   \sum _{ \bb{p} \in \mathcal{L}^\ast } f(\bb{p}) \phi _{-\bb{p}}( \bb{A} \bb{\omega} )
=    \sum _{ \bb{n} \in \mathbb{Z}^d } f( \bb{A}^\ast \bb{n}) \, e^{-2\pi \ii \bb{\omega}^\Tr \bb{A}^{\ast \Tr} \bb{A} \bb{n} } \\
=&    \sum _{ \bb{n} \in \mathbb{Z}^d } f(\bb{A}^\ast \bb{n})  \, e^{-2\pi \ii \bb{\omega}^\Tr \bb{n} } = \sum _{ \bb{n} \in \mathbb{Z}^d } f[\bb{n}]  \, e^{-2\pi \ii \bb{\omega}^\Tr \bb{n} },\end{aligned}$$

where we defined $f[\bb{n}] =  f( \bb{A}^\ast \bb{n})$. In what follows,
with some abuse of notation, the discrete sequence of samples,
$f| _{\mathcal{L}^\ast} = \{ f _{\bb{n}} = f[\bb{n}] \} _{ \bb{n} \in \mathbb{Z}^d }$,
as well as its particular element will be denoted simply as $f _{\bb{n}}$
or $f[\bb{n}]$. Square brackets will emphasize that we are working in
the digital (i.e., discetely sampled) domain. We will refer to
$f[\bb{n}]$ as to a *discrete signal* and to $F[\bb{\omega}]$ as to its
*(discrete-domain) Fourier transform*. By convention, $F[\bb{\omega}]$
is extended periodically to entire $\mathbb{R}^d$.

Note that, as before, the Fourier transform is defined as the orthogonal
projection,

$$F[\bb{n}] = \langle  f,  \phi _{\bb{\omega}} \rangle _{\ell^2(\mathbb{Z}^d)},$$

only this time the inner product is on the space of sequences

$$\langle  f,  g \rangle _{\ell^2(\mathbb{Z}^d)} = \sum _{\bb{n} \in \mathbb{Z}^d} f[\bb{n}] g^\ast[\bb{n}].$$

The inverse transform is given by $$\begin{aligned}
f[\bb{n}] = \int _{\left[ -\frac{1}{2}, \frac{1}{2} \right]^d} F[\bb{\omega}] \, e^{2\pi \ii \bb{\omega}^\Tr \bb{n} } d\bb{\omega}.\end{aligned}$$
Since the discrete-domain Fourier transform is defined through its
(periodized) continuous counterpart, all the properties of the latter
are straightforwardly inherited.

## Discrete systems


Let $f$ be a continuous-domain signal, $h$ the impulse response of an
LSI system, and let us denote $g = f \ast h$. Sampling $g$ on the
lattice $\mathcal{L}^\ast = \bb{A}^\ast \mathbb{Z}^d$ yields

$$\bed _{\mathcal{L}^\ast}  \cdot (f \ast h) \fff  \mathrm{vol}\, \mathcal{L} \cdot \bed _{\mathcal{L}} \ast (F \cdot H)$$

If either $f$ or $h$ is $\mathcal{L} _0$-band limited, so is $f \ast h$
and we can further write
$\bed _{\mathcal{L}} \ast (F \cdot H) = (\bed _{\mathcal{L}} \ast F) \cdot (\bed _{\mathcal{L}} \ast H)$,
hence 

$$\begin{aligned}
\bed _{\mathcal{L}^\ast}  \cdot g =& \bed _{\mathcal{L}^\ast}  \cdot (f \ast h) = (\bed _{\mathcal{L}^\ast} \cdot f) \ast  (\bed _{\mathcal{L}^\ast} \cdot h) = \left( \sum _{\bb{p} \in \mathcal{L}^\ast} f(\bb{p}) \tau _{\bb{p}} \delta \right) \ast 
 \left( \sum _{\bb{p} \in \mathcal{L}^\ast} h(\bb{p}) \tau _{\bb{p}} \delta \right) \\
 =& \sum _{\bb{p}, \bb{q} \in \mathcal{L}^\ast} f(\bb{p}) h(\bb{q})  \tau _{\bb{p}+\bb{q}} \delta 
 = \sum _{\bb{s} \in \mathcal{L}^\ast} \left( \sum _{\bb{p} \in \mathcal{L}^\ast}  f(\bb{p}) h(\bb{s}-\bb{p})  \right) \tau _{\bb{s}} \delta \\
 =& \sum _{\bb{m} \in \mathbb{Z}^d} \left( \sum _{\bb{n} \in \mathbb{Z}^d }  f(\bb{A}^\ast \bb{n}) h(\bb{A}^\ast(\bb{m}-\bb{n}))  \right) \tau _{\bb{A}^\ast \bb{m}} \delta \\
 =&  \sum _{\bb{m} \in \mathbb{Z}^d} \left( \sum _{\bb{n} \in \mathbb{Z}^d }  f[ \bb{n}] h[\bb{m}-\bb{n}]  \right) \tau _{\bb{A}^\ast \bb{m}} \delta.\end{aligned}$$

On the other hand,

$$\bed _{\mathcal{L}^\ast}  \cdot g = \sum _{\bb{m} \in \mathbb{Z}^d  } g[\bb{m}] \tau _{\bb{A}^\ast \bb{m}} \delta.$$

Hence, we arrive at a discrete form of convolution

$$g[\bb{n}] = \sum _{\bb{m} \in \mathbb{Z}^d }  f[ \bb{m}] h[\bb{n}-\bb{m}] = (f \ast h)[\bb{n}].$$

It is straightforward to show that convolution is transformed to
point-wise product in the frequency domain.

### Kroenecker’s delta

Recall that Dirac’s delta was a distribution that acted like identity
under convolution. The discrete-domain analog is the *Kroenecker delta*,
which is the sequence $\delta[\bb{n}]$ defined as $\delta[\bb{0}]=1$ and
$\delta[\bb{n}]=0$ for every $\bb{n} \ne \bb{0}$.

## Decimation


Let $f$ be a signal sampled on the lattice
$\mathcal{L}^\ast = \bb{A} \mathbb{Z}^d$. From the definition of the
Fourier transform, the following relation holds between the discrete and
continuous domains: 

$$\begin{aligned}
F[\bb{n}] = \sum _{ \bb{n} \in \mathbb{Z}^d } f(\bb{A}^\ast  \bb{n})  \, e^{-2\pi \ii \bb{\omega}^\Tr \bb{n} } 
=  \mathrm{vol}\, \mathcal{L} \, \sum _{\bb{n} \in \mathbb{Z}^d } F(\bb{A}   (\bb{\omega} - \bb{n})).\end{aligned}$$

Let $\bb{M}^\ast$ be an integer $d \times d$ matrix such that
$\mathcal{M}^\ast = \bb{M}^\ast \mathbb{Z}^d \subset \mathbb{Z}^d$ (that
is, the lattice $\mathcal{M}^\ast$ is has only integer points). We
denote by $\mathcal{M}^\ast _0$ the (discrete) *structural element* of
$\mathcal{M}^\ast$, i.e., the smallest set such that

$$\bigcup _{\bb{p} \in \mathcal{M}^\ast } (\bb{p} + \mathcal{M}^\ast _0)  = \mathbb{Z}^d.$$

Let us *sub-sample* $f[\bb{n}]$ on $\mathcal{M}^\ast$ creating a new
discrete signal
$g[\bb{n}] = f[\bb{M}^\ast\bb{n}] = f(\bb{A}^\ast \bb{M}^\ast \bb{n})$
(this operation is ofter referred to as *compression*). Essentially,
$g[\bb{n}]$ can be thought of $f(\bb{x})$ sampled on the dilated lattice
$\mathcal{L}^{\prime \ast} = \bb{A}^\ast \bb{M}^\ast \mathbb{Z}^d $
(with the corresponding dual
$\mathcal{L}^{\prime} = \bb{A} \bb{M} \mathbb{Z}^d$, where
$\bb{M}^\ast = \bb{M}^{-\Tr}$). Using definitions and elementary
algebra, we obtain the following expression for the Fourier transform of
$g[\bb{n}]$: 

$$\begin{aligned}
G[\bb{\omega}] =& \sum _{ \bb{n} \in \mathbb{Z}^d } f(\bb{A}^\ast \bb{M}^\ast \bb{n})  \, e^{-2\pi \ii \bb{\omega}^\Tr \bb{n} } 
=  \mathrm{vol}\, \mathcal{L}^{\prime} \, \sum _{\bb{n} \in \mathbb{Z}^d } F(\bb{A}  \bb{M}  (\bb{\omega} - \bb{n}))  \\
=&  \mathrm{vol}\, \mathcal{L} \cdot  \mathrm{vol}\, \mathcal{M}  \, \sum _{\bb{n} \in \mathbb{Z}^d } \sum _{\bb{m} \in \mathcal{M}^\ast _0 } F(\bb{A}  ( \bb{M} (\bb{\omega} - \bb{m} ) -\bb{n} ) ) \\
=& \mathrm{vol}\, \mathcal{M}    \, \sum _{\bb{m} \in  \mathcal{M}^\ast _0 }  \mathrm{vol}\, \mathcal{L} \, \sum _{\bb{n} \in \mathbb{Z}^d }    F(\bb{A}  ( \bb{M}  (\bb{\omega} - \bb{m} ) -\bb{n} ) ) \\
=&  \mathrm{vol}\, \mathcal{M}   \, \sum _{\bb{m} \in \mathcal{M}^\ast _0 }  F[\bb{M} (\bb{\omega} - \bb{m}) ].\end{aligned}$$


Firstly, notice that the continous-domain Fourier transform of $f$ does
not appear in the final expression; $G[\bb{\omega}]$ is entirely
expressed in terms of the discrete-domain Fourier transform
$F[\bb{\omega}]$. Therefore, we do not need to assume any underlying
continuos-domain signal $f(\bb{x})$ (though it can always be obtained at
least in principle via band-limited interpolation of $f[\bb{n}]$), and
perform the sub-sampling of $f[\bb{n}]$ entirely in the discrete domain.

Secondly, observe the striking resemblance of the above expression to
its continuous counterpart. The discrete-domain Fourier transform of the
signal being sampled is periodized on the dual lattice $\mathcal{M}$
(recall that sub-sampling is performed on $\mathcal{M}^\ast$). As in the
continuous case, the support of $F[\bb{\omega}]$ has to be contained
within the structural element of the dual lattice $\mathcal{M}$ in order
to avoid aliasing; in other words, $f[\bb{n}]$ has to be
$\mathcal{M}$-band limited. In that case, the replicas of
$F[\bb{\omega}]$ will not overlap, and $f[\bb{n}]$ can be recovered
exactly from $g[\bb{n}]$ by keeping only the main replica of
$F[\bb{\omega}]$ in $G[\bb{\omega}]$ and suppressing the others.

### Anti-aliasing

In general, when the signal is not band-limited and an anti-aliasing
filter has to be applied *before* sub-sampling. The frequency response
of the ideal low-pass filter $H _{\mathcal{M}}[\bb{\omega}]$ has to
satisfy $\mathrm{supp}(H _{\mathcal{M}}) \subset \mathcal{M} _0$ (of
course, this is a slightly sloppy notation: $H _{\mathcal{M}}$ is
periodic, hence the support is replicated on $\mathcal{M}$). Such a
low-pass filter is achieved by
$H _{\mathcal{M}}[\bb{\omega}] = \mathrm{box}( \bb{M}^{-1} \bb{\omega} )$
(again, periodized on $\mathcal{M}$), whose impulse response is
straightforwardly obtained by invoking the stretching property of the
Fourier transform

$$h _{\mathcal{M}}[\bb{n}] =  \mathrm{vol}\,\mathcal{M} \cdot \mathrm{sinc}(\bb{M}^\Tr \bb{n}).$$

The combined action of digital anti-aliasing filtering followed by
sub-sampling is usually referred to as *decimation*. On block diagrams,
we will denote sub-sampling as $\downarrow _{\bb{M}^\ast}$ or
$\bb{M}^\ast : \bb{I}$, followed by a decimation filter
$H _{\mathcal{M}}$.

## Interpolation


Let $f$ be a signal sampled, as before, on $\mathcal{L}^\ast $, and let
$\bb{N}^\ast$ be a $d \times d$ matrix such that
$\mathbb{Z}^d \subset \mathcal{N}^\ast = \bb{N}^\ast \mathbb{Z}^d $
(that is, a sub-sample of the lattice $\mathcal{N}^\ast$ is the integer
lattice). We define the *expansion* of $f[\bb{n}]$ to $\mathcal{N}^\ast$
as

$$g[\bb{n}] = \left\{ \begin{array}{ll} f[\bb{N}^\ast \bb{n}] & : \bb{N}^\ast \bb{n} \in \bb{Z}^d  \\
0 & : \mathrm{otherwise}.
\end{array}\right.$$ 

Substituing $\bb{p} = \bb{N}^\ast \bb{n}$ and,
inversely, $\bb{n} = \bb{N}^\Tr \bb{p}$, we obtain

$$g[\bb{n}] = \sum _{ \bb{p} \in \mathbb{Z}^d } f[\bb{p}] \delta[\bb{n} - \bb{N}^\Tr \bb{p}].$$

(observe the resemblance to the continuous case). In the frequency
domain, this translates to 

$$\begin{aligned}
G[\bb{\omega}] =& \sum _{ \bb{n} \in \mathbb{Z}^d } g[\bb{n}]  \, e^{-2\pi \ii \bb{\omega}^\Tr \bb{n} } = \sum _{ \bb{p} \in \mathbb{Z}^d } f[\bb{p}]  \, e^{-2\pi \ii \bb{\omega}^\Tr \bb{N}^\Tr \bb{p} } 
= F[\bb{N}\bb{\omega}].\end{aligned}$$ 

Since $\mathrm{vol}\, \mathcal{N}^\ast < 1$, the periodic discrete-time
Fourier transform $F[\bb{\omega}]$ is shrunk, and a few of its periods
alias to higher frequencies of $G[\bb{\omega}]$. However, unlike the
case of sub-sampling where information was lost, all the information
about $f[\bb{n}]$ is still present in $g[\bb{n}]$. In order to fully
recover it, we just need to remove the superfluous replicas of
$F[\bb{\omega}]$ outside
$\bb{N}^{-1} \left[ -\frac{1}{2}, \frac{1}{2} \right]^d$. The frequency
response of the ideal low-pass filter $H _{\mathcal{N}}[\bb{\omega}]$ has
to satisfy
$\mathrm{supp}(H _{\mathcal{N}}) \subset \bb{N}^{-1} \left[ -\frac{1}{2}, \frac{1}{2} \right]^d$.

Denoting $\mathcal{L}^{\prime} = \bb{A} \bb{N} \mathbb{Z}^d$, by
definition of the discrete-domain Fourier transform, 

$$\begin{aligned}
F[\bb{N}\bb{\omega}] H _{\mathcal{N}}[\bb{\omega}]
=& \mathrm{vol}\, \mathcal{L} \, \sum _{\bb{n} \in \mathbb{Z}^d } F(\bb{A}   (\bb{N} \bb{\omega} - \bb{n})) H _{\mathcal{N}}[\bb{\omega}] \\
=& \frac{\mathrm{vol}\, \mathcal{L}'}{\mathrm{vol}\, \mathcal{N}} \, \sum _{\bb{n} \in \mathbb{Z}^d } F(\bb{A} \bb{N}  ( \bb{\omega} - \bb{N}^{-1} \bb{n}))  H _{\mathcal{N}}[\bb{\omega}] \\
=&  \mathrm{vol}\, \mathcal{L}' \, \sum _{\bb{p} \in \mathbb{Z}^d } F(\bb{A} \bb{N}  ( \bb{\omega} - \bb{p})),
\end{aligned}
$$ 

since $H _{\mathcal{N}}(\bb{\omega})$ removes all replicas of $F$
except those with integer $\bb{N}^{-1}\bb{n}$. Note that we also assumed
that the filter has DC gain of $\mathrm{vol}\, \mathcal{N}$ to make the
result consistent with the sampling of $f(\bb{x})$ on
$\mathcal{L}^{\prime \ast}$. Such a filter is achieved by
$H _{\mathcal{N}}(\bb{\omega}) =  \mathrm{vol}\, \mathcal{N} \, \mathrm{box}( \bb{N} \bb{\omega} )$.
The corresponding impulse response is

$$h _{\mathcal{N}}[\bb{n}] =  \mathrm{sinc}(\bb{N}^\ast \bb{n}).$$

Convolving $h _{\mathcal{N}}[\bb{n}]$ with $g[\bb{n}]$ yields the
following identity 

$$\begin{aligned}
f[\bb{n}] =& (g \ast h _{\mathcal{N}})[\bb{n}]  = \sum _{ \bb{p} \in \mathbb{Z}^d } f[\bb{p}] \, h _{\mathcal{N}}[\bb{n} - \bb{N}^\Tr \bb{p}]  =  \sum _{ \bb{m} \in \mathbb{Z}^d } f[\bb{N}^\ast \bb{m}] \, \mathrm{sinc}(\bb{N}^\ast \bb{n} -  \bb{m}).\end{aligned}$$

The combined action of expansion followed by anti-aliasing filtering is
called *interpolation*. Note the striking resemblance to the continuous
case. As in the continuous case, here as well practice mandates the use
of other interpolation kernels that have finite support, which in the
frequency domain is translated to non-ideal low-pass filters. On block
diagrams, we will denote expansion by $\uparrow _{\bb{N}^\ast}$ or
$\bb{I} : \bb{N}^\ast$, followed by an interpolation filter
$H _{\mathcal{N}}$.

## Resampling


Let $f$ be a signal sampled, as before, on $\mathcal{L}^\ast $, and let
$\mathcal{Q}^\ast = \bb{Q}^\ast \mathbb{Z}^d$ be a new lattice onto
which we would like to *resample* $f$ without leaving the digital
domain. If $\bb{Q}^\ast$ can be factorized into
$\bb{Q}^\ast = \bb{M}^\ast \bb{N}^\ast$ such that
$\mathbb{Z}^d \subset  \bb{N}^\ast \mathbb{Z}^d$ and
$\bb{M}^\ast \mathbb{Z}^d \subset \mathbb{Z}^d$, we can first
interpolate $f[\bb{n}]$ to $\bb{N}^\ast \mathbb{Z}^d$ and then decimate
it onto $\bb{M}^\ast \mathbb{Z}^d$ as described above. On block diagrams
we will denote resamplers as $\bb{M}^\ast : \bb{N}^\ast$. They are
particularly useful when resampling images to a new grid related to the
original one by a non-integer (yet, rational) factor.

A few elementary rules can be straightforwardly proved:

1.  *Successive sub-sampling cascades can be combined*:
    $\downarrow _{\bb{M}^\ast _1} \longrightarrow\, \downarrow _{\bb{M}^\ast _2} \, =\,  \downarrow _{\bb{M} _2^\ast \bb{M} _1^\ast}$

2.  *Successive expansion cascades can be combined*:
    $\uparrow _{\bb{N}^\ast _1} \longrightarrow\, \uparrow _{\bb{N}^\ast _2} \, =\, \uparrow _{\bb{N} _2^\ast \bb{N} _1^\ast}$

3.  *Expansion can be exactly inverted*:
    $\uparrow _{\bb{N}^\ast} \longrightarrow\, \downarrow _{\bb{N}^{-\ast}} \, =\, \mathrm{id}$

4.  *Sub-sampling cannot usually be exactly inverted*:
    $\downarrow _{\bb{M}^\ast} \longrightarrow\, \uparrow _{\bb{M}^{-\ast}} \, \ne \, \mathrm{id}$


## Noble identities

Treating sub-sampling and expansion and discrete-domain systems, let us
quickly highlight several of their properties.

### Linearity

Observe that the sub-sampling operation we defined,
$(\downarrow _{\bb{M}^\ast} f)[\bb{n}] = f[\bb{M}^\ast \bb{n}]$ is
linear, i.e., it commutes with addition and scaling:

$$(\downarrow _{\bb{M}^\ast} (\alpha f + \beta g))[\bb{n}]  = (\alpha f + \beta g) [\bb{M}^\ast \bb{n}] = \alpha f [\bb{M}^\ast \bb{n}] + \beta g [\bb{M}^\ast \bb{n}] = \alpha (\downarrow _{\bb{M}^\ast} f)[\bb{n}] + \beta (\downarrow _{\bb{M}^\ast} g)[\bb{n}].$$

The same is also true for the expansion operation

$$(\uparrow _{\bb{N}^\ast} f)[\bb{n}] = \left\{ \begin{array}{ll} f[\bb{N}^\ast \bb{n}] & : \bb{N}^\ast \bb{n} \in \bb{Z}^d  \\
0 & : \mathrm{otherwise},
\end{array}\right.$$ 

since 

$$\begin{aligned}
(\uparrow _{\bb{N}^\ast} (\alpha f + \beta g))[\bb{n}] =&  \left\{ \begin{array}{ll} (\alpha f + \beta g)[\bb{N}^\ast \bb{n}] & : \bb{N}^\ast \bb{n} \in \bb{Z}^d  \\
0 & : \mathrm{otherwise}
\end{array}\right.  \\
=& 
\alpha \left\{ \begin{array}{ll} f [\bb{N}^\ast \bb{n}] & : \bb{N}^\ast \bb{n} \in \bb{Z}^d  \\
0 & : \mathrm{otherwise}
\end{array}\right. + \beta \left\{ \begin{array}{ll} g [\bb{N}^\ast \bb{n}] & : \bb{N}^\ast \bb{n} \in \bb{Z}^d  \\
0 & : \mathrm{otherwise}.
\end{array}\right. \\
=&
\alpha (\uparrow _{\bb{N}^\ast} f)[\bb{n}] + \beta (\uparrow _{\bb{N}^\ast} g)[\bb{n}].\end{aligned}$$


### Translation

When changing the order of translation and sub-sampling, the lattice
matrix multiplies the translation vector:

$$(\tau _{\bb{p}} \downarrow _{\bb{M}^\ast}  f)[\bb{n}] =   f[\bb{M}^\ast (\bb{n} - \bb{p})] =   f[\bb{M}^\ast \bb{n} - \bb{M}^\ast \bb{p}] = ( \downarrow _{\bb{M}^\ast} \tau _{\bb{M}^\ast \bb{p}}  f)[\bb{n}].$$

Similarly, for the expansion operation, 

$$\begin{aligned}
( \uparrow _{\bb{N}^\ast} \tau _{\bb{p}}  f)[\bb{n}] =&   \left\{ \begin{array}{ll} f[\bb{N}^\ast \bb{n}-\bb{p}] & : \bb{N}^\ast \bb{n}   \in \bb{Z}^d  \\
0 & : \mathrm{otherwise},
\end{array}\right.
=
 \left\{ \begin{array}{ll} f[\bb{N}^\ast (\bb{n}-\bb{N}^{-\ast} \bb{p})] & : \bb{N}^\ast \bb{n}   \in \bb{Z}^d  \\
0 & : \mathrm{otherwise},
\end{array}\right. \\
=&  \left\{ \begin{array}{ll} f[\bb{N}^\ast (\bb{n}-\bb{N}^{-\ast} \bb{p})] & : \bb{N}^\ast (\bb{n} - \bb{N}^{-\ast} \bb{p})   \in \bb{Z}^d  \\
0 & : \mathrm{otherwise},
\end{array}\right. \\
=& ( \uparrow _{\bb{N}^\ast}   f)[\bb{n} -  \bb{N}^{-\ast} \bb{p}] = (\tau _{ \bb{N}^{-\ast} \bb{p}}  \uparrow _{\bb{N}^\ast}  f)[\bb{n}].\end{aligned}$$

Note that since the translation operator is modified when swapped with
the sub-sampling/expansion, these systems are linear but *not* shift
invariant!

### Exchanging sub-sampling and filtering

Let $\mathcal{H}$ be an LSI system defined by 

$$\begin{aligned}
(\mathcal{H}f)[\bb{n}] =& (h \ast f)[\bb{n}] = \sum _{\bb{m} \in \mathbb{Z}^d} h[\bb{m}] f[\bb{n}-\bb{m}]  = \left( \sum _{\bb{m} \in \mathbb{Z}^d} h[\bb{m}] \delta[\bb{n}-\bb{m}] \right) \ast f[\bb{n}] \\
=& \left( \sum _{\bb{m} \in \mathbb{Z}^d} h[\bb{m}] \tau _{\bb{m}} \right) f[\bb{n}]\end{aligned}$$

Combining linearity and translation properties of sub-sampling, we
obtain 

$$\begin{aligned}
 \mathcal{H} \downarrow _{\bb{M}^\ast}  =& \sum _{\bb{m} \in \mathbb{Z}^d} h[\bb{m}] ( \tau _{\bb{m}}  \downarrow _{\bb{M}^\ast} )
= \sum _{\bb{m} \in \mathbb{Z}^d} h[\bb{m}]  (\downarrow _{\bb{M}^\ast} \tau _{\bb{M}^\ast \bb{m}} )  \\
=& \downarrow _{\bb{M}^\ast}  \left( \sum _{\bb{m} \in \mathbb{Z}^d} h[\bb{m}] \tau _{\bb{M}^\ast \bb{m}}   \right) = \downarrow _{\bb{M}^\ast} \mathcal{G}.\end{aligned}$$

In other words, applying the system $\mathcal{H}$ followed by
sub-sampling on $\mathcal{M}^\ast$ is equivalent to first sub-sampling
on $\mathcal{M}^\ast$ and then applying the system

$$\mathcal{G} =  \sum _{\bb{n} \in \mathbb{Z}^d} h[\bb{n}] \tau _{\bb{M}^\ast \bb{n}} = \sum _{\bb{p} \in \bb{M}^\ast \mathbb{Z}^d} h[\bb{M}^{-\ast}\bb{p}] \tau _{\bb{p}},$$

which can be described by the *expanded* impulse response

$$g[\bb{n}] = \left\{ \begin{array}{ll} h [\bb{M}^{-\ast} \bb{n}] & : \bb{M}^{-\ast} \bb{n} \in \bb{Z}^d  \\
0 & : \mathrm{otherwise} 
\end{array}\right.$$ 

The frequency response of this equivalent system is
straightforwardly given by 

$$\begin{aligned}
G[\bb{\omega}] =&   \sum _{\bb{n} \in \mathbb{Z}^d} h[\bb{n}]\, e^{-2\pi \ii \, \bb{\omega}^\Tr \bb{M}^\ast \bb{n}} =  \sum _{\bb{n} \in \mathbb{Z}^d} h[\bb{n}]\, e^{-2\pi \ii \, (\bb{M}^{\ast \Tr} \bb{\omega})^\Tr \bb{n}} = H[\bb{M}^{\ast \Tr} \bb{\omega}] = H[\bb{M}^{-1} \bb{\omega}].\end{aligned}$$

This result is usually known as a *noble identity* and can be summarized
as 

$$\begin{aligned}
  h \ast (\downarrow _{\bb{M}^\ast} f) \, = \,  \downarrow _{\bb{M}^\ast}  (( \uparrow _{\bb{M}^{-\ast}} h  ) \ast f )   \end{aligned}$$

or, alternatively, as the following informal expression

$$\begin{aligned}
    H[\bb{\omega}] \downarrow _{\bb{M}^\ast}  \, = \, \downarrow _{\bb{M}^\ast} H[\bb{M}^{-1} \bb{\omega}]   \end{aligned}$$

where with some abuse of notation we represent the action of an LSI
system by its Fourier transform.

### Exchanging expansion and filtering

Similarly to the previous case, combining linearity and translation
properties of expansion, we obtain 

$$\begin{aligned}
 \uparrow _{\bb{N}^\ast} \mathcal{H}  =& \sum _{\bb{m} \in \mathbb{Z}^d} h[\bb{m}] (\uparrow _{\bb{N}^\ast}  \tau _{\bb{m}}   )
= \sum _{\bb{m} \in \mathbb{Z}^d} h[\bb{m}]  ( \tau _{\bb{N}^{-\ast} \bb{m}} \uparrow _{\bb{N}^\ast} )  \\
=&  \left( \sum _{\bb{m} \in \mathbb{Z}^d} h[\bb{m}] \tau _{\bb{N}^{-\ast} \bb{m}}   \right) \uparrow _{\bb{N}^\ast} =  \mathcal{G} \uparrow _{\bb{N}^\ast}.\end{aligned}$$

with

$$\mathcal{G} =  \sum _{\bb{n} \in \mathbb{Z}^d} h[\bb{n}] \tau _{\bb{N}^{-\ast} \bb{n}} = \sum _{\bb{p} \in \bb{N}^{-\ast} \mathbb{Z}^d} h[\bb{N}^{\ast}\bb{p}] \tau _{\bb{p}},$$

that can be described by the expanded impulse response

$$g[\bb{n}] = \left\{ \begin{array}{ll} h [\bb{N}^{\ast} \bb{n}] & : \bb{N}^{\ast} \bb{n} \in \bb{Z}^d  \\
0 & : \mathrm{otherwise} 
\end{array}\right.$$ 

The frequency response of this equivalent system is
straightforwardly given by 

$$\begin{aligned}
G[\bb{\omega}] =&   \sum _{\bb{n} \in \mathbb{Z}^d} h[\bb{n}]\, e^{-2\pi \ii \, \bb{\omega}^\Tr \bb{N}^{-\ast} \bb{n}} =  \sum _{\bb{n} \in \mathbb{Z}^d} h[\bb{n}]\, e^{-2\pi \ii \, (\bb{N} \bb{\omega})^\Tr \bb{n}} =  H[\bb{N} \bb{\omega}].\end{aligned}$$

This yields our second *noble identity* that can be summarized as

$$\begin{aligned}
  \uparrow _{\bb{N}^\ast} ( h \ast  f) \, = \,   (\uparrow _{\bb{N}^{\ast}} h ) \ast  (\uparrow _{\bb{N}^\ast} f)    \end{aligned}$$

or, alternatively, as the following informal expression

$$\begin{aligned}
   \uparrow _{\bb{N}^\ast}  H[\bb{\omega}]  \, = \, H[\bb{N} \bb{\omega}] \uparrow _{\bb{N}^\ast}    \end{aligned}$$

## Polyphase representation

Let us fix some lattice
$\mathcal{M}^\ast = \bb{M}^\ast \mathbb{Z}^d \subset \mathbb{Z}^d$ and
denote by $\mathcal{M}^\ast _0$ its structural element. In this notation,
we can represent

$$\sum _{\bb{n} \in \mathbb{Z}^d} f[\bb{n}] = \sum _{\bb{m} \in \mathcal{M}^\ast _0 } \sum _{\bb{n} \in \mathbb{Z}^d} f[\bb{M}^\ast \bb{n} + \bb{m}].$$

Let $\mathcal{H}$ be an LSI system with the impulse response $h[\bb{n}]$
that can be described as

$$\mathcal{H} =  \sum _{\bb{n} \in \mathbb{Z}^d} h[\bb{n}] \tau _{\bb{n}}.$$

We define the *polyphase components* of $h[\bb{n}]$ with respect to
$\mathcal{M}^\ast$ as

$$e _{\bb{m}}[\bb{n}] = h[\bb{M}^\ast \bb{n} + \bb{m}]$$ 

for every $\bb{m} \in \mathcal{M}^\ast _0$. Each such component can be regarded as
the impulse response of a system in its own right,

$$\mathcal{E} _{\bb{m}} =  \sum _{\bb{n} \in \mathbb{Z}^d} e _{\bb{m}}[\bb{n}] \tau _{\bb{n}} 
=  \sum _{\bb{n} \in \mathbb{Z}^d} h[\bb{M}^\ast \bb{n} + \bb{m}] \tau _{\bb{n}}.$$

In these terms, the system $\mathcal{H}$ itself can be expressed as

$$\begin{aligned}
\mathcal{H} =&   \sum _{\bb{m} \in \mathcal{M}^\ast _0 } \sum _{\bb{n} \in \mathbb{Z}^d}  h[\bb{M}^\ast \bb{n} + \bb{m}] \tau _{\bb{M}^\ast \bb{n} + \bb{m}} =  \sum _{\bb{m} \in \mathcal{M}^\ast } \left( \sum _{\bb{n} \in \mathbb{Z}^d}  e _{\bb{m}}[\bb{n}] \tau _{\bb{M}^\ast \bb{n}} \right) \tau _{\bb{m}}.\end{aligned}$$

Note that the inner sum is just an expanded impulse response of
$\mathcal{E} _{\bb{m}}$,

$$(h \ast f)[\bb{n}] =  \sum _{\bb{m} \in \mathcal{M}^\ast }  (\uparrow _{\bb{M}^{-\ast}} e _{\bb{m}} ) \ast f[\bb{n}-\bb{m}]$$

which in the frequency domain assumes the form of

$$H[\bb{\omega}] =  \sum _{\bb{m} \in \mathcal{M}^\ast _0 } E _{\bb{m}}[\bb{M}^{-1} \bb{\omega}]  \, e^{-2\pi \ii \, \bb{\omega}^\Tr \bb{m}}.$$

This decomposition of an LSI system in terms of the polyphase components
is known as the *type I polyphase decomposition*.

Alternatively, we can define the polyphase sequences
$r _{\bb{m}}[\bb{n}] = e _{-\bb{m}}[\bb{n}] = h[\bb{M}^\ast \bb{n} - \bb{m}]$,
yielding the *type II polyphase decomposition*:

$$\mathcal{H} =  \sum _{\bb{m} \in \mathcal{M}^\ast } \tau _{-\bb{m}}  \left( \sum _{\bb{n} \in \mathbb{Z}^d}  r _{\bb{m}}[\bb{n}] \tau _{\bb{M}^\ast \bb{n}} \right) .$$

## Efficient implementation of decimators


When $\mathcal{H}$ is a system implementing a decimation filter, it is
followed by the sub-sampling $\downarrow _{\bb{M}^\ast}$. Substituting
type I polyphase decomposition with respect to $\mathcal{M}^\ast$ and
applying the noble identity yields 

$$\begin{aligned}
\downarrow _{\bb{M}^\ast} (h \ast f) =& \downarrow _{\bb{M}^\ast}  \left( \sum _{\bb{m} \in \mathcal{M}^\ast }   (\uparrow _{\bb{M}^{-\ast}} e _{\bb{m}} ) \ast (\tau _{\bb{m}} f) \right)
=
 \sum _{\bb{m} \in \mathcal{M}^\ast }  \downarrow _{\bb{M}^\ast}  (\uparrow _{\bb{M}^{-\ast}} e _{\bb{m}} ) \ast (\tau _{\bb{m}} f)  \\
 =& 
  \sum _{\bb{m} \in \mathcal{M}^\ast }  e _{\bb{m}}  \ast ( \downarrow _{\bb{M}^\ast} \tau _{\bb{m}} f) 
  =
   \sum _{\bb{m} \in \mathcal{M}^\ast }  \mathcal{E} _{\bb{m}} \downarrow _{\bb{M}^\ast} \tau _{\bb{m}} f.\end{aligned}$$

or in operator writing 

$$\begin{aligned}
\downarrow _{\bb{M}^\ast} \mathcal{H} =\sum _{\bb{m} \in \mathcal{M}^\ast }  \mathcal{E} _{\bb{m}} \downarrow _{\bb{M}^\ast} \tau _{\bb{m}} \end{aligned}$$

Both sides of the identity can be viewed as two different
implementations of the decimator system.

Let us compare the computational complexity of these two
implementations. For simplicity, let us assume $\bb{M}^\ast = M \bb{I}$,
so that $\mathcal{M}^\ast _0 = \{0,1,\dots,M-1\}^d$, and furthermore, let
us assume that $h[\bb{n}]$ is supported on the discrete square
$\{0,\dots,K-1\}^d$, while the input signal is supported on
$\{0,\dots,N-1\}^d$. To simplify expressions, let us further assume that
both $N$ and $K$ are divisible by $M$. The first system first calculated
$N^d$ outputs of the filter $\mathcal{H}$ (containing $K^d$
coefficients), each of which requiring $K^d$ multiply-and-accumulate
(MACC) operations. This totals in $(KN)^d$ MACC operations. If the input
signals are images coming at a rate of one image per second, this
requires a processor running at $(KN)^d$ hertz. For $d=2$, $M=4$, $K=1$,
and $N=4000$ (a $16$ megapixel image decimated $\times$ in each axis
using a $4 \times 4$ bicubic kernel), this amounts to $256$ MHz! Note
that afterwards, only $1/16$-th of these output samples are retained by
the sub-sampler; the rest are thrown away.

The second implementation first filters $M^d$ shifted and sub-sampled
versions of the input signal (containing
$\displaystyle{ \left(\frac{N}{M}\right)^d }$ samples) with the filters
$\mathcal{E} _{\bb{m}}$ (containing
$\displaystyle{ \left(\frac{K}{M}\right)^d }$ coefficients). Each filter
output sample requires only
$\displaystyle{ \left(\frac{KN}{M^2}\right)^d }$ MACC operations per
filter, and
$\displaystyle{ M^d \left(\frac{KN}{M^2}\right)^d = \left(\frac{KN}{M}\right)^d}$
MACC operations in total. Note the save up by a factor of $M^d$. For the
same settings as before, we would need a processor running only at $16$
MHz. Also note that the $M^d$ filters can be computed in parallel on
$M^d$ cores running at $M^{-2d}$ the frequency of the core required to
implement the first system. In our example, if the processor had $16$
cores, we could run each core at $1$ MHz. Since processor power
consumption is roughly proportional to $f^3$ (or worse), the second
system requires roughly $M^{5d}$ less power. In our example, $16$
$1$-MHz cores would consume about one million times less power. For this
reason, it is not an exaggeration to say that *any* ASIC implementing
signal decimation uses the polyphase implementation.

## Efficient implementation of interpolators

We will now repeat the same analysis for the case when $\mathcal{H}$ is
an interpolation filter preceding expansion $\uparrow _{\bb{N}^\ast}$.
Substituting the type II polyphase representation with respect to
$\mathcal{N}^{-\ast}$ and applying the noble identity yields

$$\begin{aligned}
 h \ast ( \uparrow _{\bb{N}^\ast} f) =&   \sum _{\bb{m} \in \mathcal{N}^{-\ast} }  \tau _{-\bb{m}}  (\uparrow _{\bb{N}^{\ast}} r _{\bb{m}} ) \ast ( \uparrow _{\bb{N}^\ast} f)   
=
 \sum _{\bb{m} \in \mathcal{N}^{-\ast} } \tau _{-\bb{m}}  \uparrow _{\bb{N}^{\ast}} ( r _{\bb{m}} \ast f)  \\
 =& 
 \sum _{\bb{m} \in \mathcal{N}^{-\ast} } \tau _{-\bb{m}}  \uparrow _{\bb{N}^{\ast}} \mathcal{R} _{\bb{m}}  f .\end{aligned}$$

or in operator writing 

$$\begin{aligned}
\mathcal{H}  \uparrow _{\bb{N}^\ast} =\sum _{\bb{m} \in \mathcal{N}^{-\ast} }  \tau _{-\bb{m}}  \uparrow _{\bb{N}^{\ast}} \mathcal{R} _{\bb{m}}   \end{aligned}$$


As before, both sides of the identity can be viewed as two different
implementations of the interpolator system. The second system first
filters the inputs with the smaller filters $\mathcal{R} _{\bb{m}}$, the
expands the results and sums their shifted versions.
