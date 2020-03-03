# dsstar
Simple demo code for the paper 

DS*: Tighter Lifting-Free Convex Relaxations for Quadratic Matching Problems
F. Bernard, C. Theobalt, M. Moeller
IEEE Conference on Computer Vision and Pattern Recognition (CVPR). June, 2018

If you use this code you are required to cite this paper.

The code provides the implementation of a convex relaxation of the quadratic assignment problem (QAP), which reads
min_X X(:)'*W*X(:) s.t. X is a permutation matrix

Here, X is an n-by-n permutation matrix, and W is an n^2-by-n^2 matrix that defines the QAP. For further details see the above paper.

