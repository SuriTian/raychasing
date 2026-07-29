# Introduction
Knots are actually really cool!!! If you had a string, and scrambled it up a bit then tied the ends, the chances are that it can be transformed into a very cool equation. This is called the "Conway Polynomial". 

This repo is a combination of ray tracing (physics to get a really pretty image) + mathematical knots. 

We take an input of a "braid" representation of a knot (e.g. [1,1,1] for trefoil) and uses skein relations to compute its conway polynomial. 

Unfortunately there is an incredibly large amount of type of knots in this world, so not all of them will be rendered. So the program checks if the knot you gave it is renderable per current standards, and renders it if it is. 

Warning: rendering a knot may take a VERY long time. This is because calculating the ppm just takes a long time, especially if we want to make it look cool. 
Calculating the conway polynomial though, won't take that long. Can refer a prototype image below to what it is supposed to look like. 
This project is still in the works to support more types of knots and improve user experience. **new patch lives in a different branch from main**

## Braid Notation

Braids are represented as a list of signed integers corresponding to the generators of the Artin braid group.

- A positive integer `k` represents the generator $\sigma_k$.
- A negative integer `-k` represents the inverse generator $\sigma_k^{-1}$.

For example:

| Braid word | Vector representation |
|------------|-----------------------|
| $\sigma_1$ | `[1]` |
| $\sigma_1^{-1}$ | `[-1]` |
| $\sigma_2$ | `[2]` |
| $\sigma_2^{-1}$ | `[-2]` |
| $\sigma_1\sigma_2^{-1}\sigma_1$ | `[1, -2, 1]` |
| $(\sigma_1\sigma_2)^3$ | `[1, 2, 1, 2, 1, 2]` |

For more information, you can read this: https://en.wikipedia.org/wiki/Braid_group

# Raw Compilation 
cmake --build build  
./build/Debug/RayChasing.exe > image.ppm  

Cannot be viewed without a ppm converter.

# KNOT PROTOTYPE
![Trefoil render](images/trefoil_knot_prototype.png)

![Command Line Calculations](images/command_line_prototype.png)