---
permalink: /info/
title: Course Info
classes: wide
---

The purpose of this course is a self-contained introduction to modern techniques in image processing. We will understand concepts involved in acquiring, sampling, representing, compressing, and processing of multi-dimensional signals (images and videos). We will discuss the notion of inverse problems and ways to solve them -- from traditional prior-based methods to modern learning-based methods. Finally, apart from the standard image acquisition model, we will also look at more exotic computational and medical imaging schemes and understand the challenges and applications involved there within.

## Learning Outcomes

At the end of the course, the student will:

1.	Understand the evolution of ideas in image processing: 
1.	Know how to implement these leading image processing algorithms in Python.
1.	Perform a small research project using the studied notions and techniques.


## Administration

**Evaluation**: 60% Homework assignments, 40% final project. (tentative)

**Language**: The course will be taught in English.

**Credits**: 3.0.

## Course Staff

{% include course_staff.html %}

<!-- ## Literature

{% include literature.html %} -->

## Detailed Syllabus

| #    | Date         | Lecture                                                                                                                                                                               | Tutorial                                                                            |
| ---- | -----------  | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -----------------------------------------------------------------------------------   |
| 1    | `22/10/2018` | Multi-dimensional signals and systems                             | Signal processing (convolution, translation) with Python 
| 2    | `28/10/2018` | Sampling and interpolation                                        | FFT, convolution theorem, aliasing effects (e.g. in MRI)   |  |
| 3    | `04/11/2018` | Discrete-domain systems                                           | Spatial & frequency domain filtering |
| 4    | `11/11/2018` | Random signals                                                    | Auto-correlation, CCF, Wiener filter                           |
| 5    | `18/11/2018` | Inverse problems, ML & MAP estimators                             | ML vs. MAP on image denoising                                  |
| 6    | `25/11/2018` | Patch-based priors                                                | Image denoising using NLM and BM3D                             |
| 7    | `09/12/2018` | L1-L2 optimization (ISTA)                                         | L1-L2 implementation                                           |
|      |              | Final project: list of papers released                            |                                                                | 
| 8    | `16/12/2018` | Shift-invariant dictionaries, LISTA                               | LISTA using autograd                                           |
| 9    | `23/12/2018` | From parsimonious models to deep learning                         | Practical introduction to deep learning                        |
| 10   | `30/12/2018` | CNNs for image restoration                                        | Image denoising and super-resolution with CNNs                 |
| 11   | `06/01/2019` | Compressed sensing, Johnson-Lindenstrauss lemma                   | Perceptual & adversarial loss for image restoration, style-transfer  |
| 12   | `13/01/2019` | Computational and medical imaging                                 | MRI and Ultrasound imaging                                                                |
| 13   | `20/01/2019` | Project presentations                                             | ----|
| ---- | -----------  | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -----------------------------------------------------------------------------------   |

