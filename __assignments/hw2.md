---
title: "Homework 2"
permalink: assignments/hw2
excerpt: Blind super resolution
author: sanketh
date: 2022-1-12
published: true
---

**Submission date**: February 22nd, 2022

## Topics

- Resampling
- Inverse filtering
- ML, MAP
- Patch-based priors
- Super resolution

## Downloading

The assignment is available
[here](https://github.com/vistalab-technion/cs236860-hw/tree/master/hw2).

We recommend you use `git` to clone the repo:
```shell
git clone https://github.com/vistalab-technion/cs236860-hw.git
```
This will allow you to pull updates from us in the event that they are needed.

## FAQ

**Q**: How can I express 2D convolutions using matrix-vector multiplication?

**A**: You can stack the columns of the image in a single vector and use a circulant matrix (with some missing rows) as your operator. See the following [post](https://stackoverflow.com/questions/16798888/2-d-convolution-as-a-matrix-matrix-multiplication).

**Q**: How can I express downsampling using matrix-vector multiplication?

**A**: You can stack the columns of the image in a single vector. Then take the identity matrix and discard the rows that correspond to the entries you want to discard from the vector. Downsampling is equivalent to the multiplication between the matrix and the vector.