//    Sequence Manipulation Suite. A collection of simple JavaScript programs
//    for generating, formatting, and analyzing short DNA and protein
//    sequences.
//    Copyright (C) 2020 Paul Stothard stothard@ualberta.ca
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <https://www.gnu.org/licenses/>.
//

// Written by Paul Stothard, University of Alberta, Canada
// Replaced pairwise alignment core with a typed-array affine-gap implementation.

(function (global) {
  "use strict";

  var NEG_INF = -1073741824;
  var SMALL_PROBLEM_SIZE = 4096;

  var STATE_M = 0;
  var STATE_X = 1;
  var STATE_Y = 2;
  var MODE_START = 3;
  var MODE_ANY = 4;

  function safeAdd(score, delta) {
    if (score <= NEG_INF / 2) {
      return NEG_INF;
    }
    return score + delta;
  }

  function safeCombine(left, right) {
    if (left <= NEG_INF / 2 || right <= NEG_INF / 2) {
      return NEG_INF;
    }
    return left + right;
  }

  function buildMatrix(name, alphabet, rows, aliases) {
    var size = alphabet.length;
    var lookup = new Int16Array(128);
    var i;

    for (i = 0; i < lookup.length; i++) {
      lookup[i] = -1;
    }

    for (i = 0; i < size; i++) {
      lookup[alphabet.charCodeAt(i)] = i;
    }

    if (aliases) {
      for (i = 0; i < aliases.length; i++) {
        lookup[aliases[i][0].charCodeAt(0)] = lookup[aliases[i][1].charCodeAt(0)];
      }
    }

    return {
      name: name,
      alphabet: alphabet,
      size: size,
      scores: new Int16Array(rows),
      lookup: lookup,
    };
  }

  function encodeSequence(sequence, scoringMatrix) {
    var upper = sequence.toUpperCase();
    var length = upper.length;
    var encoded = new Uint8Array(length);
    var display = new Uint16Array(length);
    var i;
    var code;
    var mapped;

    for (i = 0; i < length; i++) {
      code = upper.charCodeAt(i);
      mapped = code < 128 ? scoringMatrix.lookup[code] : -1;
      if (mapped < 0) {
        throw (
          "Unsupported residue '" +
          upper.charAt(i) +
          "' for matrix " +
          scoringMatrix.name +
          "."
        );
      }
      encoded[i] = mapped;
      display[i] = code;
    }

    return {
      encoded: encoded,
      display: display,
    };
  }

  function initializeRow(targetM, targetX, targetY, startMode, gapOpen, gapExtend) {
    var j;

    for (j = 0; j < targetM.length; j++) {
      targetM[j] = NEG_INF;
      targetX[j] = NEG_INF;
      targetY[j] = NEG_INF;
    }

    if (startMode === MODE_START || startMode === STATE_M) {
      targetM[0] = 0;
    } else if (startMode === STATE_X) {
      targetX[0] = 0;
    } else if (startMode === STATE_Y) {
      targetY[0] = 0;
    } else if (startMode === MODE_ANY) {
      targetM[0] = 0;
      targetX[0] = 0;
      targetY[0] = 0;
    } else {
      throw "Unknown alignment state.";
    }

    for (j = 1; j < targetM.length; j++) {
      targetY[j] = Math.max(
        safeAdd(targetM[j - 1], -gapOpen),
        safeAdd(targetY[j - 1], -gapExtend)
      );
    }
  }

  function computeForwardRowScores(seqOne, seqTwo, startMode, scoringMatrix, gapOpen, gapExtend) {
    var n = seqTwo.length;
    var prevM = new Int32Array(n + 1);
    var prevX = new Int32Array(n + 1);
    var prevY = new Int32Array(n + 1);
    var currM = new Int32Array(n + 1);
    var currX = new Int32Array(n + 1);
    var currY = new Int32Array(n + 1);
    var i;
    var j;
    var rowOffset;
    var residueOne;
    var tmp;

    initializeRow(prevM, prevX, prevY, startMode, gapOpen, gapExtend);

    for (i = 1; i <= seqOne.length; i++) {
      for (j = 0; j <= n; j++) {
        currM[j] = NEG_INF;
        currX[j] = NEG_INF;
        currY[j] = NEG_INF;
      }

      currX[0] = Math.max(
        safeAdd(prevM[0], -gapOpen),
        safeAdd(prevX[0], -gapExtend)
      );

      residueOne = seqOne[i - 1];
      rowOffset = residueOne * scoringMatrix.size;

      for (j = 1; j <= n; j++) {
        currX[j] = Math.max(
          safeAdd(prevM[j], -gapOpen),
          safeAdd(prevX[j], -gapExtend)
        );
        currY[j] = Math.max(
          safeAdd(currM[j - 1], -gapOpen),
          safeAdd(currY[j - 1], -gapExtend)
        );

        tmp = prevM[j - 1];
        if (prevX[j - 1] > tmp) {
          tmp = prevX[j - 1];
        }
        if (prevY[j - 1] > tmp) {
          tmp = prevY[j - 1];
        }

        currM[j] = safeAdd(tmp, scoringMatrix.scores[rowOffset + seqTwo[j - 1]]);
      }

      tmp = prevM;
      prevM = currM;
      currM = tmp;

      tmp = prevX;
      prevX = currX;
      currX = tmp;

      tmp = prevY;
      prevY = currY;
      currY = tmp;
    }

    return {
      m: prevM,
      x: prevX,
      y: prevY,
    };
  }

  function computeReverseRowScores(seqOne, seqTwo, endMode, scoringMatrix, gapOpen, gapExtend) {
    var n = seqTwo.length;
    var prevM = new Int32Array(n + 1);
    var prevX = new Int32Array(n + 1);
    var prevY = new Int32Array(n + 1);
    var currM = new Int32Array(n + 1);
    var currX = new Int32Array(n + 1);
    var currY = new Int32Array(n + 1);
    var i;
    var j;
    var rowOffset;
    var residueOne;
    var tmp;

    initializeRow(prevM, prevX, prevY, endMode, gapOpen, gapExtend);

    for (i = 1; i <= seqOne.length; i++) {
      for (j = 0; j <= n; j++) {
        currM[j] = NEG_INF;
        currX[j] = NEG_INF;
        currY[j] = NEG_INF;
      }

      currX[0] = Math.max(
        safeAdd(prevM[0], -gapOpen),
        safeAdd(prevX[0], -gapExtend)
      );

      residueOne = seqOne[seqOne.length - i];
      rowOffset = residueOne * scoringMatrix.size;

      for (j = 1; j <= n; j++) {
        currX[j] = Math.max(
          safeAdd(prevM[j], -gapOpen),
          safeAdd(prevX[j], -gapExtend)
        );
        currY[j] = Math.max(
          safeAdd(currM[j - 1], -gapOpen),
          safeAdd(currY[j - 1], -gapExtend)
        );

        tmp = prevM[j - 1];
        if (prevX[j - 1] > tmp) {
          tmp = prevX[j - 1];
        }
        if (prevY[j - 1] > tmp) {
          tmp = prevY[j - 1];
        }

        currM[j] = safeAdd(
          tmp,
          scoringMatrix.scores[rowOffset + seqTwo[seqTwo.length - j]]
        );
      }

      tmp = prevM;
      prevM = currM;
      currM = tmp;

      tmp = prevX;
      prevX = currX;
      currX = tmp;

      tmp = prevY;
      prevY = currY;
      currY = tmp;
    }

    return {
      m: prevM,
      x: prevX,
      y: prevY,
    };
  }

  function chooseBestEndState(scoreM, scoreX, scoreY, endMode) {
    if (endMode === STATE_M) {
      return {
        state: STATE_M,
        score: scoreM,
      };
    }
    if (endMode === STATE_X) {
      return {
        state: STATE_X,
        score: scoreX,
      };
    }
    if (endMode === STATE_Y) {
      return {
        state: STATE_Y,
        score: scoreY,
      };
    }

    if (scoreX > scoreM) {
      if (scoreY > scoreX) {
        return {
          state: STATE_Y,
          score: scoreY,
        };
      }
      return {
        state: STATE_X,
        score: scoreX,
      };
    }

    if (scoreY > scoreM) {
      return {
        state: STATE_Y,
        score: scoreY,
      };
    }

    return {
      state: STATE_M,
      score: scoreM,
    };
  }

  function solveSmallAffine(seqOne, displayOne, seqTwo, displayTwo, startMode, endMode, scoringMatrix, gapOpen, gapExtend) {
    var rows = seqOne.length + 1;
    var cols = seqTwo.length + 1;
    var total = rows * cols;
    var scoresM = new Int32Array(total);
    var scoresX = new Int32Array(total);
    var scoresY = new Int32Array(total);
    var traceM = new Uint8Array(total);
    var traceX = new Uint8Array(total);
    var traceY = new Uint8Array(total);
    var i;
    var j;
    var idx;
    var diagIndex;
    var upIndex;
    var leftIndex;
    var scoreFromM;
    var scoreFromX;
    var scoreFromY;
    var best;
    var rowOffset;
    var state;
    var alignedOne = [];
    var alignedTwo = [];

    for (idx = 0; idx < total; idx++) {
      scoresM[idx] = NEG_INF;
      scoresX[idx] = NEG_INF;
      scoresY[idx] = NEG_INF;
      traceM[idx] = 255;
      traceX[idx] = 255;
      traceY[idx] = 255;
    }

    if (startMode === MODE_START || startMode === STATE_M) {
      scoresM[0] = 0;
    } else if (startMode === STATE_X) {
      scoresX[0] = 0;
    } else if (startMode === STATE_Y) {
      scoresY[0] = 0;
    } else if (startMode === MODE_ANY) {
      scoresM[0] = 0;
      scoresX[0] = 0;
      scoresY[0] = 0;
    }

    for (j = 1; j < cols; j++) {
      idx = j;
      leftIndex = j - 1;
      scoreFromM = safeAdd(scoresM[leftIndex], -gapOpen);
      scoreFromY = safeAdd(scoresY[leftIndex], -gapExtend);

      if (scoreFromM >= scoreFromY) {
        scoresY[idx] = scoreFromM;
        traceY[idx] = STATE_M;
      } else {
        scoresY[idx] = scoreFromY;
        traceY[idx] = STATE_Y;
      }
    }

    for (i = 1; i < rows; i++) {
      idx = i * cols;
      upIndex = idx - cols;
      scoreFromM = safeAdd(scoresM[upIndex], -gapOpen);
      scoreFromX = safeAdd(scoresX[upIndex], -gapExtend);

      if (scoreFromM >= scoreFromX) {
        scoresX[idx] = scoreFromM;
        traceX[idx] = STATE_M;
      } else {
        scoresX[idx] = scoreFromX;
        traceX[idx] = STATE_X;
      }

      rowOffset = seqOne[i - 1] * scoringMatrix.size;

      for (j = 1; j < cols; j++) {
        idx = i * cols + j;
        diagIndex = idx - cols - 1;
        upIndex = idx - cols;
        leftIndex = idx - 1;

        scoreFromM = safeAdd(scoresM[upIndex], -gapOpen);
        scoreFromX = safeAdd(scoresX[upIndex], -gapExtend);
        if (scoreFromM >= scoreFromX) {
          scoresX[idx] = scoreFromM;
          traceX[idx] = STATE_M;
        } else {
          scoresX[idx] = scoreFromX;
          traceX[idx] = STATE_X;
        }

        scoreFromM = safeAdd(scoresM[leftIndex], -gapOpen);
        scoreFromY = safeAdd(scoresY[leftIndex], -gapExtend);
        if (scoreFromM >= scoreFromY) {
          scoresY[idx] = scoreFromM;
          traceY[idx] = STATE_M;
        } else {
          scoresY[idx] = scoreFromY;
          traceY[idx] = STATE_Y;
        }

        best = scoresM[diagIndex];
        state = STATE_M;
        if (scoresX[diagIndex] > best) {
          best = scoresX[diagIndex];
          state = STATE_X;
        }
        if (scoresY[diagIndex] > best) {
          best = scoresY[diagIndex];
          state = STATE_Y;
        }

        scoresM[idx] = safeAdd(best, scoringMatrix.scores[rowOffset + seqTwo[j - 1]]);
        traceM[idx] = state;
      }
    }

    idx = total - 1;
    best = chooseBestEndState(scoresM[idx], scoresX[idx], scoresY[idx], endMode);
    state = best.state;
    i = seqOne.length;
    j = seqTwo.length;

    while (i > 0 || j > 0) {
      idx = i * cols + j;

      if (state === STATE_M) {
        alignedOne.push(String.fromCharCode(displayOne[i - 1]));
        alignedTwo.push(String.fromCharCode(displayTwo[j - 1]));
        state = traceM[idx];
        i--;
        j--;
      } else if (state === STATE_X) {
        alignedOne.push(String.fromCharCode(displayOne[i - 1]));
        alignedTwo.push("-");
        state = traceX[idx];
        i--;
      } else {
        alignedOne.push("-");
        alignedTwo.push(String.fromCharCode(displayTwo[j - 1]));
        state = traceY[idx];
        j--;
      }
    }

    alignedOne.reverse();
    alignedTwo.reverse();

    return {
      alignedOne: alignedOne.join(""),
      alignedTwo: alignedTwo.join(""),
      score: best.score,
    };
  }

  function solveAffineRecursive(seqOne, displayOne, seqTwo, displayTwo, startMode, endMode, scoringMatrix, gapOpen, gapExtend, chunksOne, chunksTwo) {
    var m = seqOne.length;
    var n = seqTwo.length;
    var baseResult;
    var midpoint;
    var forward;
    var reverse;
    var j;
    var scoreM;
    var scoreX;
    var scoreY;
    var bestState = STATE_M;
    var bestSplit = 0;
    var bestScore = NEG_INF;
    var leftScore;
    var rightScore;

    if (m === 0 || n === 0 || m * n <= SMALL_PROBLEM_SIZE) {
      baseResult = solveSmallAffine(
        seqOne,
        displayOne,
        seqTwo,
        displayTwo,
        startMode,
        endMode,
        scoringMatrix,
        gapOpen,
        gapExtend
      );
      chunksOne.push(baseResult.alignedOne);
      chunksTwo.push(baseResult.alignedTwo);
      return baseResult.score;
    }

    midpoint = m >>> 1;
    forward = computeForwardRowScores(
      seqOne.subarray(0, midpoint),
      seqTwo,
      startMode,
      scoringMatrix,
      gapOpen,
      gapExtend
    );
    reverse = computeReverseRowScores(
      seqOne.subarray(midpoint),
      seqTwo,
      endMode,
      scoringMatrix,
      gapOpen,
      gapExtend
    );

    for (j = 0; j <= n; j++) {
      scoreM = safeCombine(forward.m[j], reverse.m[n - j]);
      if (scoreM > bestScore) {
        bestScore = scoreM;
        bestState = STATE_M;
        bestSplit = j;
      }

      scoreX = safeCombine(forward.x[j], reverse.x[n - j]);
      if (scoreX > bestScore) {
        bestScore = scoreX;
        bestState = STATE_X;
        bestSplit = j;
      }

      scoreY = safeCombine(forward.y[j], reverse.y[n - j]);
      if (scoreY > bestScore) {
        bestScore = scoreY;
        bestState = STATE_Y;
        bestSplit = j;
      }
    }

    leftScore = solveAffineRecursive(
      seqOne.subarray(0, midpoint),
      displayOne.subarray(0, midpoint),
      seqTwo.subarray(0, bestSplit),
      displayTwo.subarray(0, bestSplit),
      startMode,
      bestState,
      scoringMatrix,
      gapOpen,
      gapExtend,
      chunksOne,
      chunksTwo
    );
    rightScore = solveAffineRecursive(
      seqOne.subarray(midpoint),
      displayOne.subarray(midpoint),
      seqTwo.subarray(bestSplit),
      displayTwo.subarray(bestSplit),
      bestState,
      endMode,
      scoringMatrix,
      gapOpen,
      gapExtend,
      chunksOne,
      chunksTwo
    );

    return leftScore + rightScore;
  }

  function alignPairwiseAffine(sequenceOne, sequenceTwo, matrixName, gapOpen, gapExtend) {
    var scoringMatrix = MATRIX_BY_NAME[matrixName];
    var encodedOne;
    var encodedTwo;
    var chunksOne = [];
    var chunksTwo = [];
    var score;

    if (!scoringMatrix) {
      throw "Unknown scoring matrix: " + matrixName + ".";
    }

    gapOpen = parseInt(gapOpen, 10);
    gapExtend = parseInt(gapExtend, 10);

    if (isNaN(gapOpen) || isNaN(gapExtend) || gapOpen < 0 || gapExtend < 0) {
      throw "Gap opening and extension penalties must be non-negative integers.";
    }

    encodedOne = encodeSequence(sequenceOne, scoringMatrix);
    encodedTwo = encodeSequence(sequenceTwo, scoringMatrix);

    score = solveAffineRecursive(
      encodedOne.encoded,
      encodedOne.display,
      encodedTwo.encoded,
      encodedTwo.display,
      MODE_START,
      MODE_ANY,
      scoringMatrix,
      gapOpen,
      gapExtend,
      chunksOne,
      chunksTwo
    );

    return {
      alignedM: chunksOne.join(""),
      alignedN: chunksTwo.join(""),
      score: score,
    };
  }

  var EDNAFULL = buildMatrix(
    "EDNAFULL",
    "ATGCSWRYKMBVHDN",
    [
      5, -4, -4, -4, -4, 1, 1, -4, -4, 1, -4, -1, -1, -1, -2,
      -4, 5, -4, -4, -4, 1, -4, 1, 1, -4, -1, -4, -1, -1, -2,
      -4, -4, 5, -4, 1, -4, 1, -4, 1, -4, -1, -1, -4, -1, -2,
      -4, -4, -4, 5, 1, -4, -4, 1, -4, 1, -1, -1, -1, -4, -2,
      -4, -4, 1, 1, -1, -4, -2, -2, -2, -2, -1, -1, -3, -3, -1,
      1, 1, -4, -4, -4, -1, -2, -2, -2, -2, -3, -3, -1, -1, -1,
      1, -4, 1, -4, -2, -2, -1, -4, -2, -2, -3, -1, -3, -1, -1,
      -4, 1, -4, 1, -2, -2, -4, -1, -2, -2, -1, -3, -1, -3, -1,
      -4, 1, 1, -4, -2, -2, -2, -2, -1, -4, -1, -3, -3, -1, -1,
      1, -4, -4, 1, -2, -2, -2, -2, -4, -1, -3, -1, -1, -3, -1,
      -4, -1, -1, -1, -1, -3, -3, -1, -1, -3, -1, -2, -2, -2, -1,
      -1, -4, -1, -1, -1, -3, -1, -3, -3, -1, -2, -1, -2, -2, -1,
      -1, -1, -4, -1, -3, -1, -3, -1, -3, -1, -2, -2, -1, -2, -1,
      -1, -1, -1, -4, -3, -1, -1, -3, -1, -3, -2, -2, -2, -1, -1,
      -2, -2, -2, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    ],
    [["U", "T"], ["X", "N"]]
  );

  var BLOSUM62 = buildMatrix(
    "BLOSUM62",
    "ARNDCQEGHILKMFPSTWYV",
    [
      4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0,
      -1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3,
      -2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3,
      -2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3,
      0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1,
      -1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2,
      -1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2,
      0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3,
      -2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3,
      -1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3,
      -1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1,
      -1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2,
      -1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1,
      -2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1,
      -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2,
      1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2,
      0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0,
      -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3,
      -2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1,
      0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4,
    ]
  );

  var BLOSUM45 = buildMatrix(
    "BLOSUM45",
    "ARNDCQEGHILKMFPSTWYV",
    [
      5, -2, -1, -2, -1, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -2, -2, 0,
      -2, 7, 0, -1, -3, 1, 0, -2, 0, -3, -2, 3, -1, -2, -2, -1, -1, -2, -1, -2,
      -1, 0, 6, 2, -2, 0, 0, 0, 1, -2, -3, 0, -2, -2, -2, 1, 0, -4, -2, -3,
      -2, -1, 2, 7, -3, 0, 2, -1, 0, -4, -3, 0, -3, -4, -1, 0, -1, -4, -2, -3,
      -1, -3, -2, -3, 12, -3, -3, -3, -3, -3, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1,
      -1, 1, 0, 0, -3, 6, 2, -2, 1, -2, -2, 1, 0, -4, -1, 0, -1, -2, -1, -3,
      -1, 0, 0, 2, -3, 2, 6, -2, 0, -3, -2, 1, -2, -3, 0, 0, -1, -3, -2, -3,
      0, -2, 0, -1, -3, -2, -2, 7, -2, -4, -3, -2, -2, -3, -2, 0, -2, -2, -3, -3,
      -2, 0, 1, 0, -3, 1, 0, -2, 10, -3, -2, -1, 0, -2, -2, -1, -2, -3, 2, -3,
      -1, -3, -2, -4, -3, -2, -3, -4, -3, 5, 2, -3, 2, 0, -2, -2, -1, -2, 0, 3,
      -1, -2, -3, -3, -2, -2, -2, -3, -2, 2, 5, -3, 2, 1, -3, -3, -1, -2, 0, 1,
      -1, 3, 0, 0, -3, 1, 1, -2, -1, -3, -3, 5, -1, -3, -1, -1, -1, -2, -1, -2,
      -1, -1, -2, -3, -2, 0, -2, -2, 0, 2, 2, -1, 6, 0, -2, -2, -1, -2, 0, 1,
      -2, -2, -2, -4, -2, -4, -3, -3, -2, 0, 1, -3, 0, 8, -3, -2, -1, 1, 3, 0,
      -1, -2, -2, -1, -4, -1, 0, -2, -2, -2, -3, -1, -2, -3, 9, -1, -1, -3, -3, -3,
      1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -3, -1, -2, -2, -1, 4, 2, -4, -2, -1,
      0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -1, -1, 2, 5, -3, -1, 0,
      -2, -2, -4, -4, -5, -2, -3, -2, -3, -2, -2, -2, -2, 1, -3, -4, -3, 15, 3, -3,
      -2, -1, -2, -2, -3, -1, -2, -3, 2, 0, 0, -1, 0, 3, -3, -2, -1, 3, 8, -1,
      0, -2, -3, -3, -1, -3, -3, -3, -3, 3, 1, -2, 1, 0, -3, -1, 0, -3, -1, 5,
    ]
  );

  var BLOSUM80 = buildMatrix(
    "BLOSUM80",
    "ARNDCQEGHILKMFPSTWYV",
    [
      5, -2, -2, -2, -1, -1, -1, 0, -2, -2, -2, -1, -1, -3, -1, 1, 0, -3, -2, 0,
      -2, 6, -1, -2, -4, 1, -1, -3, 0, -3, -3, 2, -2, -4, -2, -1, -1, -4, -3, -3,
      -2, -1, 6, 1, -3, 0, -1, -1, 0, -4, -4, 0, -3, -4, -3, 0, 0, -4, -3, -4,
      -2, -2, 1, 6, -4, -1, 1, -2, -2, -4, -5, -1, -4, -4, -2, -1, -1, -6, -4, -4,
      -1, -4, -3, -4, 9, -4, -5, -4, -4, -2, -2, -4, -2, -3, -4, -2, -1, -3, -3, -1,
      -1, 1, 0, -1, -4, 6, 2, -2, 1, -3, -3, 1, 0, -4, -2, 0, -1, -3, -2, -3,
      -1, -1, -1, 1, -5, 2, 6, -3, 0, -4, -4, 1, -2, -4, -2, 0, -1, -4, -3, -3,
      0, -3, -1, -2, -4, -2, -3, 6, -3, -5, -4, -2, -4, -4, -3, -1, -2, -4, -4, -4,
      -2, 0, 0, -2, -4, 1, 0, -3, 8, -4, -3, -1, -2, -2, -3, -1, -2, -3, 2, -4,
      -2, -3, -4, -4, -2, -3, -4, -5, -4, 5, 1, -3, 1, -1, -4, -3, -1, -3, -2, 3,
      -2, -3, -4, -5, -2, -3, -4, -4, -3, 1, 4, -3, 2, 0, -3, -3, -2, -2, -2, 1,
      -1, 2, 0, -1, -4, 1, 1, -2, -1, -3, -3, 5, -2, -4, -1, -1, -1, -4, -3, -3,
      -1, -2, -3, -4, -2, 0, -2, -4, -2, 1, 2, -2, 6, 0, -3, -2, -1, -2, -2, 1,
      -3, -4, -4, -4, -3, -4, -4, -4, -2, -1, 0, -4, 0, 6, -4, -3, -2, 0, 3, -1,
      -1, -2, -3, -2, -4, -2, -2, -3, -3, -4, -3, -1, -3, -4, 8, -1, -2, -5, -4, -3,
      1, -1, 0, -1, -2, 0, 0, -1, -1, -3, -3, -1, -2, -3, -1, 5, 1, -4, -2, -2,
      0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -2, -1, -1, -2, -2, 1, 5, -4, -2, 0,
      -3, -4, -4, -6, -3, -3, -4, -4, -3, -3, -2, -4, -2, 0, -5, -4, -4, 11, 2, -3,
      -2, -3, -3, -4, -3, -2, -3, -4, 2, -2, -2, -3, -2, 3, -4, -2, -2, 2, 7, -2,
      0, -3, -4, -4, -1, -3, -3, -4, -4, 3, 1, -3, 1, -1, -3, -2, 0, -3, -2, 4,
    ]
  );

  var PAM70 = buildMatrix(
    "PAM70",
    "ARNDCQEGHILKMFPSTWYV",
    [
      5, -4, -2, -1, -4, -2, -1, 0, -4, -2, -4, -4, -3, -6, 0, 1, 1, -9, -5, -1,
      -4, 8, -3, -6, -5, 0, -5, -6, 0, -3, -6, 2, -2, -7, -2, -1, -4, 0, -7, -5,
      -2, -3, 6, 3, -7, -1, 0, -1, 1, -3, -5, 0, -5, -6, -3, 1, 0, -6, -3, -5,
      -1, -6, 3, 6, -9, 0, 3, -1, -1, -5, -8, -2, -7, -10, -4, -1, -2, -10, -7, -5,
      -4, -5, -7, -9, 9, -9, -9, -6, -5, -4, -10, -9, -9, -8, -5, -1, -5, -11, -2, -4,
      -2, 0, -1, 0, -9, 7, 2, -4, 2, -5, -3, -1, -2, -9, -1, -3, -3, -8, -8, -4,
      -1, -5, 0, 3, -9, 2, 6, -2, -2, -4, -6, -2, -4, -9, -3, -2, -3, -11, -6, -4,
      0, -6, -1, -1, -6, -4, -2, 6, -6, -6, -7, -5, -6, -7, -3, 0, -3, -10, -9, -3,
      -4, 0, 1, -1, -5, 2, -2, -6, 8, -6, -4, -3, -6, -4, -2, -3, -4, -5, -1, -4,
      -2, -3, -3, -5, -4, -5, -4, -6, -6, 7, 1, -4, 1, 0, -5, -4, -1, -9, -4, 3,
      -4, -6, -5, -8, -10, -3, -6, -7, -4, 1, 6, -5, 2, -1, -5, -6, -4, -4, -4, 0,
      -4, 2, 0, -2, -9, -1, -2, -5, -3, -4, -5, 6, 0, -9, -4, -2, -1, -7, -7, -6,
      -3, -2, -5, -7, -9, -2, -4, -6, -6, 1, 2, 0, 10, -2, -5, -3, -2, -8, -7, 0,
      -6, -7, -6, -10, -8, -9, -9, -7, -4, 0, -1, -9, -2, 8, -7, -4, -6, -2, 4, -5,
      0, -2, -3, -4, -5, -1, -3, -3, -2, -5, -5, -4, -5, -7, 7, 0, -2, -9, -9, -3,
      1, -1, 1, -1, -1, -3, -2, 0, -3, -4, -6, -2, -3, -4, 0, 5, 2, -3, -5, -3,
      1, -4, 0, -2, -5, -3, -3, -3, -4, -1, -4, -1, -2, -6, -2, 2, 6, -8, -4, -1,
      -9, 0, -6, -10, -11, -8, -11, -10, -5, -9, -4, -7, -8, -2, -9, -3, -8, 13, -3, -10,
      -5, -7, -3, -7, -2, -8, -6, -9, -1, -4, -4, -7, -7, 4, -9, -5, -4, -3, 9, -5,
      -1, -5, -5, -5, -4, -4, -4, -3, -4, 3, 0, -6, 0, -5, -3, -3, -1, -10, -5, 6,
    ]
  );

  var PAM30 = buildMatrix(
    "PAM30",
    "ARNDCQEGHILKMFPSTWYV",
    [
      6, -7, -4, -3, -6, -4, -2, -2, -7, -5, -6, -7, -5, -8, -2, 0, -1, -13, -8, -2,
      -7, 8, -6, -10, -8, -2, -9, -9, -2, -5, -8, 0, -4, -9, -4, -3, -6, -2, -10, -8,
      -4, -6, 8, 2, -11, -3, -2, -3, 0, -5, -7, -1, -9, -9, -6, 0, -2, -8, -4, -8,
      -3, -10, 2, 8, -14, -2, 2, -3, -4, -7, -12, -4, -11, -15, -8, -4, -5, -15, -11, -8,
      -6, -8, -11, -14, 10, -14, -14, -9, -7, -6, -15, -14, -13, -13, -8, -3, -8, -15, -4, -6,
      -4, -2, -3, -2, -14, 8, 1, -7, 1, -8, -5, -3, -4, -13, -3, -5, -5, -13, -12, -7,
      -2, -9, -2, 2, -14, 1, 8, -4, -5, -5, -9, -4, -7, -14, -5, -4, -6, -17, -8, -6,
      -2, -9, -3, -3, -9, -7, -4, 6, -9, -11, -10, -7, -8, -9, -6, -2, -6, -15, -14, -5,
      -7, -2, 0, -4, -7, 1, -5, -9, 9, -9, -6, -6, -10, -6, -4, -6, -7, -7, -3, -6,
      -5, -5, -5, -7, -6, -8, -5, -11, -9, 8, -1, -6, -1, -2, -8, -7, -2, -14, -6, 2,
      -6, -8, -7, -12, -15, -5, -9, -10, -6, -1, 7, -8, 1, -3, -7, -8, -7, -6, -7, -2,
      -7, 0, -1, -4, -14, -3, -4, -7, -6, -6, -8, 7, -2, -14, -6, -4, -3, -12, -9, -9,
      -5, -4, -9, -11, -13, -4, -7, -8, -10, -1, 1, -2, 11, -4, -8, -5, -4, -13, -11, -1,
      -8, -9, -9, -15, -13, -13, -14, -9, -6, -2, -3, -14, -4, 9, -10, -6, -9, -4, 2, -8,
      -2, -4, -6, -8, -8, -3, -5, -6, -4, -8, -7, -6, -8, -10, 8, -2, -4, -14, -13, -6,
      0, -3, 0, -4, -3, -5, -4, -2, -6, -7, -8, -4, -5, -6, -2, 6, 0, -5, -7, -6,
      -1, -6, -2, -5, -8, -5, -6, -6, -7, -2, -7, -3, -4, -9, -4, 0, 7, -13, -6, -3,
      -13, -2, -8, -15, -15, -13, -17, -15, -7, -14, -6, -12, -13, -4, -14, -5, -13, 13, -5, -15,
      -8, -10, -4, -11, -4, -12, -8, -14, -3, -6, -7, -9, -11, 2, -13, -7, -6, -5, 10, -7,
      -2, -8, -8, -8, -6, -7, -6, -5, -6, 2, -2, -9, -1, -8, -6, -6, -3, -15, -7, 7,
    ]
  );

  var MATRIX_BY_NAME = {
    ednafull: EDNAFULL,
    blosum62: BLOSUM62,
    blosum45: BLOSUM45,
    blosum80: BLOSUM80,
    pam70: PAM70,
    pam30: PAM30,
  };

  global.PairwiseAffineAlignment = {
    align: alignPairwiseAffine,
  };

  if (typeof module !== "undefined" && module.exports) {
    module.exports = global.PairwiseAffineAlignment;
  }
})(this);
