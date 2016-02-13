# Global





* * *

### code(inputCodewords, inputAlphabetSize) 

**Parameters**

**inputCodewords**: `Array`, the codewords.

**inputAlphabetSize**: `Number`, the size of the alphabet.

**Returns**: `code`, the code.


## Class: code


### code.codewords() 

**Returns**: `Array`, the codewords.

### code.codeword(i) 

**Parameters**

**i**: `Number`, the index.

**Returns**: `Array`, the codeword at index `i`.

### code.alphabetSize() 

**Returns**: `Number`, the size of the alphabet.

### code.clone() 

**Returns**: `code`, a copy of the code.

### code.length() 

**Returns**: `Number`, the length of a codeword.

### code.size() 

**Returns**: `Number`, the number of codewords in the code.

### code.informationRate() 

**Returns**: `Number`, the information rate of the code.

### code.distanceFrom(codeword2, ignoreSimilar, metric) 

**Parameters**

**codeword2**: `Array`, another codeword.

**ignoreSimilar**: `Boolean`, ignores codeword in the code if it is the same as `codeword2`.

**metric**: `String`, the distance metric used. The Hamming distance is used if no string is specified.

**Returns**: `Number`, the distance of the code from `codeword2`.

### code.weight(metric) 

**Parameters**

**metric**: `String`, the distance metric used. The Hamming distance is used if no string is specified.

**Returns**: `Number`, the weight of the code.

### code.maximumInnerProductWith(code2, type) 

**Parameters**

**code2**: `code`, another code.

**type**: `String`, the type of inner product used. The Euclidean inner product is used if no string is specified.

**Returns**: `Number`, the maximum inner product of the code and `code2`.

### code.isOrthogonalTo(code2, type) 

**Parameters**

**code2**: `code`, another code.

**type**: `String`, the type of inner product used. The Euclidean inner product is used if no string is specified.

**Returns**: `Boolean`, the maximum inner product of the code and `code2` is 0.

### code.contains(code2) 

**Parameters**

**code2**: `code`, another code.

**Returns**: `Boolean`, all codewords in `code2` is in the code.

### code.isLinear() 

**Returns**: `Boolean`, the code is linear.

### code.dimension() 

**Returns**: `code`, the dimension of the code if it is linear; `undefined` otherwise.

### code.distance(metric) 

**Parameters**

**metric**: `String`, the distance metric used. The Hamming distance is used if no string is specified.

**Returns**: `Number`, the distance of the code.

### code.relativeMinimumDistance(metric) 

**Parameters**

**metric**: `String`, the distance metric used. The Hamming distance is used if no string is specified.

**Returns**: `Number`, the relative minimum distance of the code.

### code.exactlyErrorDetecting(metric) 

**Parameters**

**metric**: `String`, the distance metric used. The Hamming distance is used if no string is specified.

**Returns**: `Number`, the number at which the code is exactly error detecting.

### code.exactlyErrorCorrecting(metric) 

**Parameters**

**metric**: `String`, the distance metric used. The Hamming distance is used if no string is specified.

**Returns**: `Number`, the number at which the code is exactly error correcting.

### code.fromBasis(basis) 

**Parameters**

**basis**: `Array`, the basis.

**Returns**: `code`, all possible codewords in the span of the basis, taking `alphabetSize()` from the code.

### code.allPossibleCodewords(length) 

**Parameters**

**length**: `Number`, the length of the code.

**Returns**: `code`, all possible codewords of length `length`.

### code.basis() 

**Returns**: `Array`, *a* basis of the linear code; `undefined` if the code is non-linear.

### code.orthogonalComplement(type) 

**Parameters**

**type**: `String`, the type of inner product used. The Euclidean inner product is used if no string is specified.

**Returns**: `code`, the orthogonal complement.

### code.dual(type) 

**Parameters**

**type**: `String`, the type of inner product used. The Euclidean inner product is used if no string is specified.

**Returns**: `code`, the dual code if the code is linear; `undefined` otherwise.

### code.isSelfOrthogonal(type) 

**Parameters**

**type**: `String`, the type of inner product used. The Euclidean inner product is used if no string is specified.

**Returns**: `boolean`, the linear code is self-orthogonal; `undefined` if the code is non-linear.

### code.isSelfDual(type) 

**Parameters**

**type**: `String`, the type of inner product used. The Euclidean inner product is used if no string is specified.

**Returns**: `boolean`, the linear code is self-dual; `undefined` if the code is non-linear.

### code.encode(basis) 

**Parameters**

**basis**: `Array`, the generator matrix.

**Returns**: `code`, the encoded code.

### code.decode(basis) 

**Parameters**

**basis**: `Array`, the generator matrix.

**Returns**: `code`, the decoded code.

### code.cosets(metric) 

**Parameters**

**metric**: `String`, the distance metric used. The Hamming distance is used if no string is specified.

**Returns**: `Array`, arrays containing a pair of elements, the first holding the array of possible coset leaders and the second holding the coset code, if the code is linear; `undefined` otherwise.

### code.syndromeLookupTable(completeDecoding, metric) 

**Parameters**

**completeDecoding**: `Boolean`, complete decoding is to be used; a random coset leader is chosen. Otherwise, use incomplete decoding, which returns `undefined` if there is more than one possible coset leader.

**metric**: `String`, the distance metric used. The Hamming distance is used if no string is specified.

**Returns**: `Array`, arrays each holding a coset leader and its syndrome code.

### code.syndromeDecoding(codeword2, completeDecoding, metric) 

**Parameters**

**codeword2**: `Array`, the codeword to be corrected.

**completeDecoding**: `Boolean`, complete decoding is to be used; a random coset leader is chosen. Otherwise, use incomplete decoding, which returns `undefined` if there is more than one possible coset leader.

**metric**: `String`, the distance metric used. The Hamming distance is used if no string is specified.

**Returns**: `Array`, the corrected codeword.

### code.lowerBound(bound, metric) 

**Parameters**

**bound**: `String`, `general`, `sphereCovering`, `singleton` or `gilbertVarshamov`. The general bound is used if no string is specified.

**metric**: `String`, the distance metric used. The Hamming distance is used if no string is specified.

**Returns**: `Number`, the lower bound.

### code.upperBound(bound, metric) 

**Parameters**

**bound**: `String`, `general`, `hamming` (or `spherePacking`), `plotkin` or `griesmer`. The general bound is used if no string is specified.

**metric**: `String`, the distance metric used. The Hamming distance is used if no string is specified.

**Returns**: `number`, the upper bound.

### code.isOptimal(metric) 

**Parameters**

**metric**: `String`, the distance metric used. The Hamming distance is used if no string is specified.

**Returns**: `Boolean`, the code is optimal.

### code.isPerfect(metric) 

**Parameters**

**metric**: `String`, the distance metric used. The Hamming distance is used if no string is specified.

**Returns**: `Boolean`, the code is perfect.

### code.isMaximumDistanceSeparable(metric) 

**Parameters**

**metric**: `String`, the distance metric used. The Hamming distance is used if no string is specified.

**Returns**: `Boolean`, the code is maximum distance separable.



* * *










