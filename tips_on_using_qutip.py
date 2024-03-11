Good Morning,

I've also had to work with outside data that are arrays/matrices. Here is what I've found to work (for what it's worth).

Qobj to numpy array / matrix:
data = q.full()
array = numpy.array(data)
matrix = numpy.matrix(data)

Numpy to Qobj
q = Qobj(array)
q = Qobj(matrix)

Numpy operations:
Matrix multiplication in numpy is '*" for matricies and (for 2d) numpy.dot(). (e.g. A.dot(B) for arrays A * B for matrices)
Element wise squaring is numpy.square() (works regardless of the dims)

Tensor on numpy things (not Qobjs):
We've been using numpy.kron(a,b) to tensor a and b (if I remember correctly this works with arrays and matrices)

Don't know if this is helpful not. :)