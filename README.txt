`fem_toolbox` is a collection of Python functions designed to perform finite element analysis (FEA) on structural systems. The package currently supports only **1D structural elements**:

- Bar elements  
- Euler-Bernoulli beam elements  
- 2D frame (beam) elements

The choice of writing an implementation from scratch instead of using well-known FEA solvers or packages comes from the fact that by designing and controlling every aspect of the code I hope to have ensured full flexibility for integrating future modules. The primary motivation for this is to **build a shape optimization module on top of this package**.

---------------------------------------------

INSTALLATION, GETTING STARTED etc: pick up the docs and start reading
