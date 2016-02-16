var code;

/**
 * @param {Array} inputCodewords the codewords represented by arrays of equal length. Each codeword must be an array of integers between *0* and *`inputAlphabetSize` - 1*, both inclusive.
 * @param {Number} inputAlphabetSize the prime size of the alphabet.
 * @returns {code} the code.
 */
code = function (inputCodewords, inputAlphabetSize) {

// Return object with key "error" if errors thrown
try{

  var range, zero, standardBasis, add, multiply, hammingDistance, euclideanInnerProduct, volume, generalLowerBound, sphereCoveringLowerBound, singletonLowerBound, generalUpperBound, hammingUpperBound, plotkinUpperBound, gilbertVarshamovLowerBound, griesmerUpperBound, inputCodewords, inputAlphabetSize, i, j, codewords, codeword, alphabetSize, clone, length, size, informationRate, distanceFrom, weight, maximumInnerProductWith, isOrthogonalTo, contains, isLinear, dimension, distance, relativeMinimumDistance, exactlyErrorDetecting, exactlyErrorCorrecting, fromBasis, allPossibleCodewords, basis, orthogonalComplement, dual, isSelfOrthogonal, isSelfDual, encode, decode, cosets, syndromeLookupTable, syndromeDecoding, lowerBound, upperBound, isOptimal, isPerfect, isMaximumDistanceSeparable;

  //// PRIVATE FUNCTIONS

  zero = function (n) {
    var outputArray, i;
    outputArray = [];
    for (i = 0; i < n; i += 1) {
      outputArray[i] = 0;
    }
    return outputArray;
  };

  standardBasis = function (length) {
    var basis, i;
    basis = [];
    for (i = 0; i < length; i += 1) {
      basis[i] = zero(length);
      basis[i][i] = 1;
    }
    return basis;
  };

  add = function (codeword, codeword2, alphabetSize) {
    var outputCodeword, i;
    outputCodeword = zero(codeword.length);
    for (i = 0; i < codeword.length; i += 1) {
      outputCodeword[i] = (codeword[i] + codeword2[i]) % alphabetSize;
    }
    return outputCodeword;
  };

  multiply = function (codeword, symbol, alphabetSize) {
    var outputCodeword, i;
    outputCodeword = zero(codeword.length);
    for (i = 0; i < codeword.length; i += 1) {
      outputCodeword[i] = (codeword[i] * symbol) % alphabetSize;
    }
    return outputCodeword;
  };

  hammingDistance = function (codeword, codeword2) {
    var d, i;
    d = 0;
    for (i = 0; i < codeword.length; i += 1) {
      if (codeword[i] !== codeword2[i]) {
        d += 1;
      }
    }
    return d;
  };

  euclideanInnerProduct = function (codeword, codeword2, alphabetSize) {
    var sum, i;
    sum = 0;
    for (i = 0; i < codeword.length; i += 1) {
      sum += (codeword[i] * codeword2[i]) % alphabetSize;
    }
    return sum % alphabetSize;
  }

  volume = function (length, alphabetSize, radius) {
    var v, i, numberOfChoices, j;
    if (radius > length) {
      return Math.pow(alphabetSize, length);
    } else {
      v = 0;
      for (i = 0; i <= radius; i += 1) {
        numberOfChoices = 1;
        for (j = 1; j <= i; j += 1) {
          numberOfChoices = numberOfChoices * (length - j + 1) / j;
        }
        v += numberOfChoices * Math.pow(alphabetSize - 1, i);
      }
      return v;
    }
  };

  generalLowerBound = function () {
    return 1;
  };

  sphereCoveringLowerBound = function (metric) {
    return Math.ceil(Math.pow(alphabetSize(), length()) / volume(length(), alphabetSize(), distance(metric) - 1));
  };

  singletonLowerBound = function (metric) {
    return Math.floor(Math.pow(alphabetSize(), (length() - distance(metric) + 1)));
  };

  generalUpperBound = function () {
    return Math.pow(alphabetSize(), length());
  };

  hammingUpperBound = function (metric) {
    return Math.floor(Math.pow(alphabetSize(), length()) / volume(length(), alphabetSize(), Math.floor((distance(metric) - 1) / 2)));
  };

  plotkinUpperBound = function (metric) {
    if (alphabetSize() === 2) {
      if (distance(metric) % 2 !== 0) {
        if (length() === 2 * distance(metric) + 1) {
          return 4 * distance(metric) + 4;
        }
        if (length() < 2 * distance(metric) + 1) {
          return 2 * Math.floor((distance(metric) + 1) / (2 * distance(metric) + 1 - length()));
        }
      }
      if (length() === 2 * distance(metric)) {
        return 4 * distance(metric);
      }
      if (length() < 2 * distance(metric)) {
        return 2 * Math.floor(distance(metric) / (2 * distance(metric) - length()));
      }
    }
    if (distance(metric) <= (1 - 1 / alphabetSize()) * length()) {
      return;
    }
    return Math.floor(distance(metric) / (distance(metric) - (1 - 1 / alphabetSize()) * length()));
  };

  gilbertVarshamovLowerBound = function (metric) {
    if (distance() < 2) {
      return;
    }
    if (volume(length() - 1, alphabetSize(), distance(metric) - 1) >= Math.pow(alphabetSize(), length() - dimension())) {
      return;
    }
    return Math.ceil(Math.pow(alphabetSize(), length() - 1) / volume(length() - 1, alphabetSize(), distance(metric) - 2));
  };

  griesmerUpperBound = function (metric) {
    var i, l;
    i = 0;
    l = Math.ceil(distance(metric) / Math.pow(alphabetSize(), i));
    while(l <= length()) {
      i += 1;
      l += Math.ceil(distance(metric) / Math.pow(alphabetSize(), i));
    }
    return Math.pow(alphabetSize(), i);
  };

  //// MODIFY THE INPUT ARGUMENTS ////

  if (!Array.isArray(inputCodewords)) {
    inputCodewords = [[]];
  }
  for (i = 0; i < inputCodewords.length; i += 1) {
    // Scalar codeword is treated as codeword of length 1
    if (typeof inputCodewords[i] === "number") {
      inputCodewords[i] = [inputCodewords[i]];
    }
    // Default codeword is empty
    if (!Array.isArray(inputCodewords[i])) {
      inputCodewords[i] = [];
    }
    // All codeword must be of the same length
    if (inputCodewords[i].length !== inputCodewords[0].length) {
      throw "Lengths of codewords differ";
    }
  }

  // Default size of alphabet is 2
  if (inputAlphabetSize !== parseInt(inputAlphabetSize, 10) || inputAlphabetSize <= 1 || !isFinite(inputAlphabetSize) || typeof inputAlphabetSize !== "number") {
    inputAlphabetSize = 2;
  }
  // Size of alphabet must be prime
  for (i = 2; i <= Math.sqrt(inputAlphabetSize); i += 1) {
    if (inputAlphabetSize % i === 0 && i !== inputAlphabetSize) {
      throw "Size of alphabet is not prime";
    }
  }

  for (i = 0; i < inputCodewords.length; i += 1) {
    // Codewords must be array of integers between 0 and alphabetSize() - 1, both inclusive
    for (j = 0; j < inputCodewords[i].length; j += 1) {
      if (inputCodewords[i][j] !== parseInt(inputCodewords[i][j], 10)){
        throw "Some symbols are not integers";
      }
      if (inputCodewords[i][j] < 0 || inputCodewords[i][j] > inputAlphabetSize - 1) {
        throw "Some symbols are not in the alphabet";
      }
    }
  }

  //// NO MORE CHANGES ARE TO BE MADE TO inputCodeword FROM HERE ON ////

  /**
   * @class code
   */

  /**
   * @returns {Array} the codewords.
   */
  codewords = function () {
    return inputCodewords.slice(0);
  };

  /**
   * @param {Number} i the index.
   * @returns {Array} the codeword at index `i`.
   */
  codeword = function (i) {
    return codewords()[i];
  };

  /**
   * @returns {Number} the size of the alphabet.
   */
  alphabetSize = function () {
    return inputAlphabetSize;
  };

  /**
   * @returns {code} a copy of the code.
   */
  clone = function () {
    return code(codewords(), alphabetSize());
  };

  /**
   * @returns {Number} the length of a codeword.
   */
  length = function () {
    return codewords()[0].length;
  };

  /**
   * @returns {Number} the number of codewords in the code.
   */
  size = function () {
    return codewords().length;
  };

  /**
   * @returns {Number} the information rate of the code.
   */
  informationRate = function () {
    return Math.log(size()) / Math.log(alphabetSize()) / length();
  };

  /**
   * @param {Array} codeword2 another codeword.
   * @param {Boolean} ignoreSimilar ignores codeword in the code if it is the same as `codeword2`.
   * @param {String} metric the distance metric used. The Hamming distance is used if no string is specified.
   * @returns {Number} the distance of the code from `codeword2`.
   */
  distanceFrom = function (codeword2, ignoreSimilar, metric) {
    var metric, d, existsSimilar, i, dist;
    if (typeof metric !== "string") {
      metric = "hamming"
    }
    d = Infinity;
    existsSimilar = false;
    if (metric === "hamming") {
      for (i = 0; i < size(); i += 1) {
        dist = hammingDistance(codeword(i), codeword2);
        if (dist === 0) {
          existsSimilar = true;
        } else {
          d = Math.min(d, dist);
        }
      }
    }
    return (ignoreSimilar || !existsSimilar) ? d : 0;
  };

  /**
   * @param {String} metric the distance metric used. The Hamming distance is used if no string is specified.
   * @returns {Number} the weight of the code.
   */
  weight = function (metric) {
    return distanceFrom(zero(length()), true, metric);
  };

  /**
   * @param {code} code2 another code.
   * @param {String} type the type of inner product used. The Euclidean inner product is used if no string is specified.
   * @returns {Number} the maximum inner product of the code and `code2`.
   */
  maximumInnerProductWith = function (code2, type) {
    var type, p, i, j;
    if (typeof type !== "string") {
      type = "euclidean"
    }
    p = 0;
    if (type === "euclidean") {
      for (i = 0; i < size(); i += 1) {
        for (j = 0; j < code2.size(); j += 1) {
          p = Math.max(p, euclideanInnerProduct(codeword(i), code2.codeword(j), alphabetSize()));
        }
      }
    }
    return p;
  };

  /**
   * @param {code} code2 another code.
   * @param {String} type the type of inner product used. The Euclidean inner product is used if no string is specified.
   * @returns {Boolean} the maximum inner product of the code and `code2` is 0.
   */
  isOrthogonalTo = function (code2, type) {
    return maximumInnerProductWith(code2, type) === 0;
  };

  /**
   * @param {code} code2 another code.
   * @returns {Boolean} all codewords in `code2` is in the code.
   */
  contains = function (code2) {
    var differ, countContainment, i, j, k;
    if (length() !== code2.length()) {
      return false;
    }
    countContainment = 0;
    for (j = 0; j < code2.size(); j += 1) {
      for (i = 0; i < size(); i += 1) {
        differ = 0;
        for (k = 0; k < length(); k += 1) {
          if (codeword(i)[k] !== code2.codeword(j)[k]) {
            differ += 1;
          }
        }
        if (differ === 0) {
          countContainment += 1;
          break;
        }
      }
    }
    return countContainment === code2.size();
  };

  /**
   * @returns {Boolean} the code is linear.
   */
  isLinear = function () {
    var lambda, mu, x, y;
    for (lambda = 0; lambda < alphabetSize(); lambda += 1) {
      for (mu = 0; mu < alphabetSize(); mu += 1) {
        for (x = 0; x < size(); x += 1) {
          for (y = 0; y < size(); y += 1) {
            if (!contains(code([add(multiply(codeword(x), lambda, alphabetSize()), multiply(codeword(y), mu, alphabetSize()), alphabetSize())], alphabetSize()))) {
              return false;
            }
          }
        }
      }
    }
    return true;
  };

  /**
   * @returns {code} the dimension of the code if it is linear; `undefined` otherwise.
   */
  dimension = function () {
    if (!isLinear()) {
      return;
    }
    return Math.log(size()) / Math.log(alphabetSize());
  };

  /**
   * @param {String} metric the distance metric used. The Hamming distance is used if no string is specified.
   * @returns {Number} the distance of the code.
   */
  distance = function (metric) {
    var metric, d, i, j;
    if (typeof metric !== "string") {
      metric = "hamming"
    }
    if (isLinear()) {
      return weight(metric);
    }
    d = Infinity;
    if (metric === "hamming") {
      for (i = 0; i < size(); i += 1) {
        for (j = i + 1; j < size(); j += 1) {
          d = Math.min(d, hammingDistance(codeword(i), codeword(j)));
        }
      }
    }
    return d;
  };

  /**
   * @param {String} metric the distance metric used. The Hamming distance is used if no string is specified.
   * @returns {Number} the relative minimum distance of the code.
   */
  relativeMinimumDistance = function (metric) {
    return (distance(metric) - 1) / length();
  };

  /**
   * @param {String} metric the distance metric used. The Hamming distance is used if no string is specified.
   * @returns {Number} the number at which the code is exactly error detecting.
   */
  exactlyErrorDetecting = function (metric) {
    return distance(metric) - 1;
  };

  /**
   * @param {String} metric the distance metric used. The Hamming distance is used if no string is specified.
   * @returns {Number} the number at which the code is exactly error correcting.
   */
  exactlyErrorCorrecting = function (metric) {
    return Math.floor((distance(metric) - 1) / 2);
  };

  /**
   * @param {Array} basis the basis.
   * @returns {code} all possible codewords in the span of the basis, taking `alphabetSize()` from the code.
   */
  fromBasis = function (basis) {
    var codewords, l, i, j, k;
    codewords = [zero(basis[0].length)];
    for (k = 0; k < alphabetSize(); k += 1) {
    }
    for (i = 0; i < basis.length; i += 1) {
      var l = codewords.length;
      for (j = 0; j < l; j += 1) {
        for (k = 0; k < alphabetSize(); k += 1) {
          codewords[l * k + j] = add(codewords[j], multiply(basis[i], k, alphabetSize()), alphabetSize());
        }
      }
    }
    return code(codewords, alphabetSize());
  };

  /**
   * @param {Number} length the length of the code.
   * @returns {code} all possible codewords of length `length`.
   */
  allPossibleCodewords = function (length) {
    return fromBasis(standardBasis(length));
  };

  /**
   * @returns {Array} *a* basis of the linear code; `undefined` if the code is non-linear.
   */
  basis = function () {
    var outputArray, outputCode, z, i;
    if (!isLinear()) {
      return;
    }
    outputArray = [];
    outputCode = code([[]], alphabetSize());
    z = code([zero(length())], alphabetSize());
    for (i = 0; i < size(); i += 1) {
      if (outputCode.contains(clone())) {
        return outputArray;
      }
      // codeword(i) must not be zero
      if (!outputCode.contains(code([codeword(i)], alphabetSize())) && !code([codeword(i)], alphabetSize()).contains(z)) {
        outputArray.push(codeword(i));
        outputCode = fromBasis(outputArray);
      }
    }
    return outputArray;
  };
  
  /**
   * @param {String} type the type of inner product used. The Euclidean inner product is used if no string is specified.
   * @returns {code} the orthogonal complement.
   */
  orthogonalComplement = function (type) {
    var code2, orthogonalCodewords, i;
    orthogonalCodewords = [];
    code2 = allPossibleCodewords(length());
    for (i = 0; i < code2.size(); i += 1) {
      if (isOrthogonalTo(code([code2.codeword(i)], alphabetSize()), type)) {
        orthogonalCodewords.push(code2.codeword(i));
      }
    }
    return code(orthogonalCodewords, alphabetSize());
  };

  /**
   * @param {String} type the type of inner product used. The Euclidean inner product is used if no string is specified.
   * @returns {code} the dual code if the code is linear; `undefined` otherwise.
   */
  dual = function (type) {
    if (!isLinear()) {
      return;
    }
    return orthogonalComplement(type);
  };

  /**
   * @param {String} type the type of inner product used. The Euclidean inner product is used if no string is specified.
   * @returns {boolean} the linear code is self-orthogonal; `undefined` if the code is non-linear.
   */
  isSelfOrthogonal = function (type) {
    if (!isLinear()) {
      return;
    }
    return dual(type).contains(clone());
  };

  /**
   * @param {String} type the type of inner product used. The Euclidean inner product is used if no string is specified.
   * @returns {boolean} the linear code is self-dual; `undefined` if the code is non-linear.
   */
  isSelfDual = function (type) {
    if (!isLinear()) {
      return;
    }
    if (!isSelfOrthogonal(type)) {
      return false;
    }
    return contains(dual(type));
  };

  /**
   * @param {Array} basis the generator matrix.
   * @returns {code} the encoded code.
   */
  encode = function (basis) {
    var outputCodewords, i, c, j;
    outputCodewords = [];
    for (i = 0; i < size(); i += 1) {
      c = zero(basis[0].length);
      for (j = 0; j < basis.length; j += 1) {
        c = add(c, multiply(basis[j], codeword(i)[j], alphabetSize()), alphabetSize());
      }
      outputCodewords.push(c);
    }
    return code(outputCodewords, alphabetSize());
  };

  /**
   * @param {Array} basis the generator matrix.
   * @returns {code} the decoded code.
   */
  decode = function (basis) {
    var inputCode, inputCodewords, k, i, c, j;
    inputCode = allPossibleCodewords(basis.length);
    inputCodewords = [];
    for (k = 0; k < size(); k += 1) {
      for (i = 0; i < inputCode.size(); i += 1) {
        c = zero(basis[0].length);
        for (j = 0; j < basis.length; j += 1) {
          c = add(c, multiply(basis[j], inputCode.codeword(i)[j], alphabetSize()), alphabetSize());
        }
        if (code([c], alphabetSize()).contains(code([codeword(k)], alphabetSize()))) {
          inputCodewords[k] = inputCode.codeword(i);
          break;
        }
      }
    }
    return code(inputCodewords, alphabetSize());
  };

  /**
   * @param {String} metric the distance metric used. The Hamming distance is used if no string is specified.
   * @returns {Array} arrays containing a pair of elements, the first holding the array of possible coset leaders and the second holding the coset code, if the code is linear; `undefined` otherwise.
   */
  cosets = function (metric) {
    var output, dist, i, C, cosetIndex, c, j, cosetCodewords, cosetCodeword, cosetLeaders;
    if (!isLinear()) {
      return;
    }
    dist = []
    for (i = 0; i < size(); i += 1) {
      dist[i] = code([codeword(i)], alphabetSize()).distanceFrom(zero(length()), false, metric);
    }
    output = [[dist, clone()]];
    C = allPossibleCodewords(length());
    cosetIndex = 1;
    for (c = 0; c < C.size(); c += 1) {
      for (j = 0; j < cosetIndex; j += 1) {
        if (output[j][1].contains(code([C.codeword(c)], alphabetSize()))) {
          break;
        }
      }
      // if C.codeword(c) cannot be found in any of the previous cosets
      if (j === cosetIndex) {
        cosetCodewords = [];
        dist = []
        for (i = 0; i < size(); i += 1) {
          cosetCodeword = add(codeword(i), C.codeword(c), alphabetSize());
          dist[i] = code([cosetCodeword], alphabetSize()).distanceFrom(zero(length()), false, metric);
          cosetCodewords[i] = cosetCodeword;
        }
        output[cosetIndex] = [dist, code(cosetCodewords, alphabetSize())];
        cosetIndex += 1;
      }
    }
    for (cosetIndex = 0; cosetIndex < output.length; cosetIndex += 1) {
      cosetLeaders = [];
      for (i = 0; i < size(); i += 1) {
        if (output[cosetIndex][0][i] === output[cosetIndex][1].distanceFrom(zero(length()), false, alphabetSize())) {
          cosetLeaders.push(output[cosetIndex][1].codeword(i));
        }
      }
      output[cosetIndex][0] = cosetLeaders;
    }
    return output;
  };

  /**
   * @param {Boolean} completeDecoding complete decoding is to be used; a random coset leader is chosen. Otherwise, use incomplete decoding, which returns `undefined` if there is more than one possible coset leader.
   * @param {String} metric the distance metric used. The Hamming distance is used if no string is specified.
   * @returns {Array} arrays each holding a coset leader and its syndrome code.
   */
  syndromeLookupTable = function (completeDecoding, metric) {
    var co, parityCheckBasis, table, cosetIndex, cosetLeader, syndrome, i, column, j;
    if (!isLinear()) {
      return;
    }
    co = cosets(metric);
    parityCheckBasis = dual().basis();
    table = [];
    for (cosetIndex = 0; cosetIndex < co.length; cosetIndex += 1) {
      if (co[cosetIndex][0].length > 1 && !completeDecoding) {
        return;
      } else {
        cosetLeader = co[cosetIndex][0][Math.floor(Math.random() * co[cosetIndex][0].length)];
        syndrome = zero(parityCheckBasis.length);
        for (i = 0; i < parityCheckBasis[0].length; i += 1) {
          column = [];
          for (j = 0; j < parityCheckBasis.length; j += 1) {
            column[j] = parityCheckBasis[j][i];
          }
          syndrome = add(syndrome, multiply(column, cosetLeader[i], alphabetSize()), alphabetSize());
        }
        table[cosetIndex] = [cosetLeader, code([syndrome], alphabetSize())];
      }
    }
    return table;
  };

  /**
   * @param {Array} codeword2 the codeword to be corrected.
   * @param {Boolean} completeDecoding complete decoding is to be used; a random coset leader is chosen. Otherwise, use incomplete decoding, which returns `undefined` if there is more than one possible coset leader.
   * @param {String} metric the distance metric used. The Hamming distance is used if no string is specified.
   * @returns {Array} the corrected codeword.
   */
  syndromeDecoding = function (codeword2, completeDecoding, metric) {
    var parityCheckBasis, syndrome, i, column, j, syndromeCode, slt, cosetLeader, additiveInverse;
    if (!isLinear()) {
      return;
    }
    parityCheckBasis = dual().basis();
    syndrome = zero(parityCheckBasis.length);
    for (i = 0; i < parityCheckBasis[0].length; i += 1) {
      column = [];
      for (j = 0; j < parityCheckBasis.length; j += 1) {
        column[j] = parityCheckBasis[j][i];
      }
      syndrome = add(syndrome, multiply(column, codeword2[i], alphabetSize()), alphabetSize());
    }
    syndromeCode = code([syndrome], alphabetSize());
    slt = syndromeLookupTable(completeDecoding, metric);
    if (typeof slt === "undefined") {
      return;
    }
    for (i = 0; i < slt.length; i += 1) {
      if (slt[i][1].contains(syndromeCode)){
        cosetLeader = slt[i][0];
        break;
      }
    }
    additiveInverse = [];
    for (i = 0; i < cosetLeader.length; i += 1) {
      additiveInverse[i] = alphabetSize() - cosetLeader[i];
    }
    return add(codeword2, additiveInverse, alphabetSize());
  };

  /**
   * @param {String} bound `general`, `sphereCovering`, `singleton` or `gilbertVarshamov`. The general bound is used if no string is specified.
   * @param {String} metric the distance metric used. The Hamming distance is used if no string is specified.
   * @returns {Number} the lower bound.
   */
  lowerBound = function (bound, metric) {
    var bound;
    if (typeof bound === "undefined") {
      bound = "general"
    }
    if (bound === "general") {
      return generalLowerBound();
    }
    if (bound === "sphereCovering") {
      return sphereCoveringLowerBound(metric);
    }
    if (bound === "singleton") {
      return singletonLowerBound(metric);
    }
    if (bound === "gilbertVarshamov") {
      if (!isLinear()) {
        return;
      }
      return gilbertVarshamovLowerBound(metric);
    }
    if (bound === "*") {
      return Math.max(lowerBound("general", metric) || 0, lowerBound("sphereCovering", metric) || 0, lowerBound("singleton", metric) || 0, lowerBound("gilberVarshamov", metric) || 0);
    }
  };

  /**
   * @param {String} bound `general`, `hamming` (or `spherePacking`), `plotkin` or `griesmer`. The general bound is used if no string is specified.
   * @param {String} metric the distance metric used. The Hamming distance is used if no string is specified.
   * @returns {number} the upper bound.
   */
  upperBound = function (bound, metric) {
    var bound;
    if (typeof bound === "undefined") {
      bound = "general"
    }
    if (bound === "general") {
      return generalUpperBound();
    }
    if (bound === "hamming" || bound === "spherePacking") {
      return hammingUpperBound(metric);
    }
    if (bound === "plotkin") {
      return plotkinUpperBound(metric);
    }
    if (bound === "griesmer") {
      if (!isLinear()) {
        return;
      }
      return griesmerUpperBound(metric);
    }
    if (bound === "*") {
      return Math.min(upperBound("general", metric) || Infinity, upperBound("hamming", metric) || Infinity, upperBound("plotkin", metric) || Infinity, upperBound("griesmer", metric) || Infinity);
    }
  };

  /**
   * @param {String} metric the distance metric used. The Hamming distance is used if no string is specified.
   * @returns {Boolean} the code is optimal.
   */
  isOptimal = function (metric) {
    if (size() === upperBound("*", metric)) {
      return true;
    }
    return undefined;
  };

  /**
   * @param {String} metric the distance metric used. The Hamming distance is used if no string is specified.
   * @returns {Boolean} the code is perfect.
   */
  isPerfect = function (metric) {
    if (size() === upperBound("hamming", metric)) {
      return true;
    }
    return undefined;
  };

  /**
   * @param {String} metric the distance metric used. The Hamming distance is used if no string is specified.
   * @returns {Boolean} the code is maximum distance separable.
   */
  isMaximumDistanceSeparable = function (metric) {
    if (!isLinear()) {
      return;
    }
    return dimension() + distance(metric) === length() + 1;
  };

  return {
    codewords: codewords,
    codeword: codeword,
    alphabetSize: alphabetSize,
    clone: clone,
    length: length,
    size: size,
    informationRate: informationRate,
    distanceFrom: distanceFrom,
    weight: weight,
    maximumInnerProductWith: maximumInnerProductWith,
    isOrthogonalTo: isOrthogonalTo,
    contains: contains,
    isLinear: isLinear,
    dimension: dimension,
    distance: distance,
    relativeMinimumDistance: relativeMinimumDistance,
    exactlyErrorDetecting: exactlyErrorDetecting,
    exactlyErrorCorrecting: exactlyErrorCorrecting,
    fromBasis: fromBasis,
    allPossibleCodewords: allPossibleCodewords,
    basis: basis,
    orthogonalComplement: orthogonalComplement,
    dual: dual,
    isSelfOrthogonal: isSelfOrthogonal,
    isSelfDual: isSelfDual,
    encode: encode,
    decode: decode,
    cosets: cosets,
    syndromeLookupTable: syndromeLookupTable,
    syndromeDecoding: syndromeDecoding,
    lowerBound: lowerBound,
    upperBound: upperBound,
    isOptimal: isOptimal,
    isPerfect: isPerfect,
    isMaximumDistanceSeparable: isMaximumDistanceSeparable,basis: basis
  };

} catch (error) {
  return {error: error};
}

};
