# polynomial.js
***a JavaScript library for working with codes (coding theory)***

License is omitted for the time being, until someone finds the files useful.

Documentation is found in the **output/** folder. It is generated using [jsdox.js](http://jsdox.org/) licensed under the MIT license.

**code.js** is the main file. Codes are created using function declaration `code(inputCodewords, inputAlphabetSize)`, not via `new`.

**tests.html** is used to check for errors using a console.

Notes:
* The symbols of the alphabet are in *0, 1, 2, ..., p*, where the size of the alphabet *p + 1* is a prime.
* Cyclic codes can be generated with the help of [polynomial.js](https://github.com/litena/polynomial.js), e.g., `code([[]],2).fromBasis(polynomial([0], "x", 2).cyclicCodeBasis(7, 3))`.

Copyright 2016 @teylim
