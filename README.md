

Efficient and stable Arnoldi restarts for matrix functions based on quadrature
================================================================================

About
-----

This MATLAB code implements a quadrature-based restarted Krylov method FUNM\_QUAD for the efficient and stable approximation of _f(A)b,_ where _A_ is a square, large and sparse matrix, _f_ is a matrix function (such that _f(A)_ is defined) and _b_ is a vector.  
  
The algorithm is described and analyzed in  

[A. Frommer, S. Güttel, and M. Schweitzer: **Efficient and stable Arnoldi restarts for matrix functions based on quadrature,** SIAM Journal on Matrix Analysis and Applications, 35(2):661--683, 2014.](https://doi.org/10.1137/13093491X)

[A. Frommer, S. Güttel, and M. Schweitzer: **Convergence of restarted Krylov subspace methods for Stieltjes functions of matrices,** SIAM Journal on Matrix Analysis and Applications, 35(4):1602--1624, 2014.](https://doi.org/10.1137/140973463)
  
[M. Schweitzer: **Restarting and error estimation in polynomial and extended Krylov subspace methods for the approximation of matrix functions,** Dissertation, Bergische Universität Wuppertal, Fakultät für Mathematik und Naturwissenschaften, 2016.](http://nbn-resolving.de/urn/resolver.pl?urn=urn:nbn:de:hbz:468-20160212-112106-7)  

 
The basic calling sequence of FUNM\_QUAD is

\[f,out\] = funm\_quad(A,b,param)

  
This MATLAB code is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. It is a research code and the authors are grateful for any constructive comments. Please refer to the above paper when using this code.

Usage
-----

To use FUNM\_QUAD extract the zip-archive on your computer and change to the directory in MATLAB. This directory also contains a PDF file with a short overview of FUNM\_QUAD, and you may also want to run the demo files to see how it works.  

Version history
---------------

*  this version: use of a relative stopping criterion in all restart cycles; controlled breakdown of Arnoldi method and automated reduction of restart length for small matrices; fixed normalisation of vector b when non-standard inner product is provided; added support for truncated (or incomplete) Arnoldi via param.truncation\_length; added a verbose parameter (0, 1, 2) to adjust level information to output on the command line; changed w = w - v\*h to w = w - v\*h(1) in Arnoldi (considerable performance gain); updated manual funm\_quad.pdf and cleaned up demo files; absorbed myintegral helper functions into a single file and removed myintegral subfolder
*  funm_quad_2015-02.zip: minor fix to prevent early breakdown in first restart cycle, relative stopping tolerance, drive\_quad.m renamed to demo\_quad.m
*  funm_quad_2014-01.zip: first version
