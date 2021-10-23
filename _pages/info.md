---
permalink: /info/
title: Course Info
classes: wide
---

The purpose of this course is a self-contained introduction to modern techniques in image processing. We will understand concepts involved in acquiring, sampling, representing, compressing, and processing of multi-dimensional signals (images and videos). We will discuss the notion of inverse problems and ways to solve them -- from traditional prior-based methods to modern learning-based methods. Finally, apart from the standard image acquisition model, we will also look at more exotic computational and medical imaging schemes and understand the challenges and applications involved there within.

## Learning Outcomes

At the end of the course, the student will:

1.  Get comfortable with the standard notions involved in acquiring, sampling, representing, compressing, and processing of images.
1.	Understand the evolution of ideas in image processing 
1.	Know how to implement these state-of-the-art image processing algorithms in Python.
1.	Perform a small research project using the studied notions and techniques.

## Administration

**Evaluation**: 60% Homework assignments, 40% final project.

**Language**: The course will be taught in English.

**Credits**: 3.0.

**Prerequisites**

This course expects **both** mathematical maturity and programming competency. The recommended prerequisites are as follows:

- A good background in linear algebra, probability and calculus. See the [supplementary material](../supplements/) page if you need a refresher on one of these.
- Programming competency. The course will be hands-on. We will use Python exclusively; it is important to have experience with it.

**Collaboration Policy and Honor Code**

By enrolling in this course, you agree that you will strictly follow our collaboration policy as specified below. Any violation of this policy will result in an immediate failure in the course, and treatment by the Technion regulations committee.

1. Submission of assignments is in singles or pairs. You are free to form study groups and discuss homeworks with other students. However, you must implement all the required code or derive all the required math independently of other groups (only with your submission partner).
2. Submitted work must be your own. You must do your own thinking, coding, debugging, and write all the answers yourself. We will run automatic plagiarism-detection software on your submissions to ensure this policy.
3. You may not use solutions from previous semesters' homeworks.
4. You may not share your solutions with other students.
5. You may not upload your homework solutions to any public website, such as Github. Private repos are OK, but they must remain so even after course completion.

## Course Staff

{% include course_staff.html %}

<!-- ## Literature

{% include literature.html %} -->

## Detailed Syllabus

The lectures will follow a flipped-classroom approach. Students will be requested to watch recorded lecture videos as a **mandatory** course requirement. We provide videos and written material on the course [Lectures](../lectures/) page to facilitate self-learning of the core topics. THe in-class lectures will be short, optional, and cover more adavanced material. In addition, any questions from the lecture videos will be addressed.

The [Tutorials](../tutorials) are based on detailed and self-contained Jupyter notebooks, which will guide you through full implementation of the concepts learned in the lecture. The in-class tutorials will cover all the material, no pre-requisite viewing is required.

The course includes hands-on homework assignments in which you'll derive and implement different image processing algorithms.

The course also includes a small research project that will require reading and presenting a research paper, implementing and improving the algorithms presented therewithin. A list of papers and further instructions will be released mid-semester.


| # | Date | Pre-requisite | Lecture | Tutorial | Homework |
| --- | --- | ---  | --- | --- | --- |
| 1 | `28/10/2021` | | Introduction | Math background | |
| 2 | `04/10/2021` | Lecture 2 | Signals & systems, Fourier | Convolutions | |
| 3 | `11/11/2021` | Lecture 3 (until Poisson summation) | Adjoint, Dirac's delta, Poisson summation | Fourier transform | |
| 4 | `18/11/2021` | Lecture 3 (the rest) | Sampling & interpolation | Computerized tomography | |
| 5 | `25/11/2021` | Lecture 4 | Discrete domain signals & systems | Sampling & iterpolation | HW1 |
| 6 | `02/11/2021` | Lecture 5a | Introduction to random signals | Inverse filtering | |
| 7 | `09/12/2021` | Lecture 5b | Inverse problems & statistical estimation | ML vs. MAP | |
| 8 | `16/12/2021` | Lecture 6 | Patch-based priors | Bilateral filters, NLM | HW2 |
| 9 | `23/12/2021` | Lecture 7 | Sparsity-based priors | PatchMatch | |
| 10 | `30/12/2021` | Lecture 8 | Structured-based priors | Sparse coding, BM3D | |
| 11 | `06/01/2022` | Lecture 9 | Learning image priors | Dictionary learning | |
| 12 | `13/01/2022` | Lecture 10 | Sampling & sensing | Deep learning (part 1) | |
| 13 | `20/01/2022` | Lecture 11 | Computational imaging | Deep learning (part 2) | HW3 |
