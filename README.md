PBF_lib
=======

c++ code for pseudo-Boolean function manipulation.

This code allows efficient read-and-write of pseudo-Boolean functions. The current implementation only includes quadratic pseudo-Boolean functions (QPBF).

We store a QPBF as a sparse matrix. This sparse matrix is stored in a heap. An coefficient $a_{ij}$ is indexed by (i,j) in the heap. 

On initialization, the user should specify the number of boolean variables in the QPBF. This is not essentially necessary, however, we make it a rule for consistency check. 

We provide two ways to access each entry $a_{ij}$. The first way is to use addTermX getTermX; the second way is via ()operator. The difference is, when an entry is accessed by ()operator, it will be created in the heap, whereas getTermX do not change the heap. Therefore, use ()operator only when the polynomial is not very large and sparsity is not necessary. 

One can completely delete an entry via delTermX, which not only set that entry to zero, but also release the memory space.

We allow +,- operation between QPBFs, as well as * operation between QPBFs and scalars. After these operations take place, some entries might become zeros, but they are still kept in the memory. One may get rid of them with clean() function. 

Clear() function destroy the whole QPBF.

size() returns the number of nonzero entries in the QPBF.

numvar() returns the number of boolean variables. 

To extend the number of boolean variables in the QPBF, one can use updatenumvar(). Currently, we do not support reducing the number of Boolean variables.

One may also evaluate the value of a QPBF with respect to a Boolean vector y using evaluate().
