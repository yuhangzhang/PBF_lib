PBF_lib
=======

c++ code for pseudo-Boolean function manipulation

This code allow efficient write and read of pseudo-Boolean fucntions. The current implementation only includes quadratic pseudo-Boolean functions (QPBF).

We store a QPBF as a sparse matrix. This sparse matrix is stored in a binary-tree, or heap. An coefficient $a_{ij}$ is indexed by (i,j) in the heap. 

On initialization, the user should specify the number of boolean variables in the QPBF. This is not essentially necessary, however, we make it a rule for consistency check. 

We provide two ways to access each entry $a_{ij}$. The first way is to use addTermX getTermX; the second way is direct access like in matlab. There is not much difference between the two, although the direct access allows editing the entries, while the other way does not.

One can completely delete an entry via delTermX, which not only set that entry to zero, but also release the memory space.

We allow +,- operation between QPBFs, as well as * operation between QPBFs and scalars. After these operations take place, some entries might become zeros, but they still kept in the memory. One way get rid of them by clean() function. 

Clear() function destroy the whole QPBF.

size() returns the number of nonzero entries in the QPBF.

numvar() returns the number of boolean variables. 

To extend or reduce (reduce is not really necessary) the number of boolean variables in the QPBF, one can use updatenumvar().

One may also evaluate the value of a QPBF with respect to a Boolean vector y using evaluate().
