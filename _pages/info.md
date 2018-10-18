---
permalink: /info/
title: Course Info
classes: wide
---

The purpose of this course is a self-contained introduction to modern techniques in image processing. We will understand concepts involved in acquiring, sampling, representing, compressing, and processing of multi-dimensional signals (images and videos). We will discuss the notion of inverse problems and ways to solve them -- from traditional prior-based methods to modern learning-based methods. Finally, apart from the standard image acquisition model, we will also look at more exotic computational and medical imaging schemes and understand the challenges and applications involved there within.

## Learning Outcomes

At the end of the course, the student will:

1.	Understand and be able to apply notions in deep learning.
1.	Know how to effectively use leading python machine-learning and deep
    learning frameworks such as PyTorch.
1.	Know how to optimize software and hardware performance in deep neural
    network applications.
1.	Know how to program GPUs using CUDA.
1.	Perform a small research project using the studied notions and techniques.


## Administration

**Evaluation**: 60% Homework assignments, 40% final project. (tentative)

**Language**: The course will be taught in English.

**Credits**: 3.0.

## Course Staff

{% include course_staff.html %}

## Literature

{% include literature.html %}

## Detailed Syllabus

| #    | Date         | Lecture                                                                                                                                                                               | Tutorial                                                                            |
| ---- | -----------  | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -----------------------------------------------------------------------------------   |
| 1    | `22/10/2018` | Multi-dimensional signals and systems                             | Signal processing (convolution, translation) with Python 
| 2    | `28/10/2018` | Sampling and interpolation                                        | FFT, convolution theorem, aliasing effects (e.g. in MRI)   |  |
| 3    | `04/11/2018` | Discrete-domain systems                                           | Spatial & frequency domain filtering |
| 4    | `11/11/2018` | Random signals                                                    | Auto-correlation, CCF, Wiener filter                                       |
| 5    | `18/11/2018` | Inverse problems, ML & MAP estimators                             | ML vs. MAP on image denoising                                             |
| 6    | `25/11/2018` | Patch-based priors                                                | Image denoising using NLM and BM3D                                                     |
| 7    | `09/12/2018` | Reinforcement learning, Deep Q-Learning, Policy Gradients (Alex)                                                                                                                      | Variational Auto Encoders                                                             |
| 8    | `16/12/2018` | CNNs on non-euclidean domains: graphs and manifolds (Alex)                                                                                                                            | Domain adaptation, style transfer                                                     |
| 9    | `23/12/2018` | Quantization and compression of DNNs (Chaim)                                                                                                                                          | Deep reinforcement learning                                                           |
| 10   | `30/12/2018` | Parallel architectures (Avi)                                                                                                                                                          | CNN on graphs                                                                         |
| 11   | `06/01/2019` | Accelerating GEMM: Custom, GPU, TPU architectures (Avi)                                                                                                                               | Introduction to CUDA 1                                                                |
| 12   | `13/01/2019` | Accelerating inference for CNNs (Avi)                                                                                                                                                 | Introduction to CUDA 2                                                                |
| 13   | `20/01/2019` | Training in distributed and parallel systems (Avi)                                                                                                                                    | CUDA kernels in PyTorch                                                               |
| ---- | -----------  | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -----------------------------------------------------------------------------------   |

