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

//Written by Paul Stothard, University of Alberta, Canada

//This class performs alignments in linear space, by recursively dividing
//the alignment. Once subalignments of acceptable size are obtained, they
//are solved using the quadratic space implementation in align_pair_quad.js.

//To use this class: (see pairwise_dna.js for example)
//var alignment = new AlignPairLinear();
//alignment.initializeMatrix(sequenceArrayM, sequenceArrayN, scoreSet);
//alignment.fillMatrix();
//alignment.align();
//var alignedSequenceStringM = alignment.getAlignedM();
//var alignedSequenceStringN = alignment.getAlignedN();

//------------------------------------ AlignPairLinear class
//AlignPairLinear class align() method
function align() {
  if (this.M.length == 0) {
    for (var j = 1; j <= this.N.length; j++) {
      this.alignedM.push("-");
      this.alignedN.push(this.N[j - 1]);
      this.score = this.score + this.scoreSet.gap;
    }
  } else if (this.N.length == 0) {
    for (var j = 1; j <= this.M.length; j++) {
      this.alignedN.push("-");
      this.alignedM.push(this.M[j - 1]);
      this.score = this.score + this.scoreSet.gap;
    }
  } else if (this.M.length == 0 && this.N.length == 0) {
    //do nothing
  } else {
    this.path(0, 0, this.M.length, this.N.length);
  }
}

//AlignPairLinear class recursive method path()
function path(i1, j1, i2, j2) {
  //alert ("i1, j1, : i2, j2 " + i1 +", " + j1 + ", " + i2 + ", " + j2);

  if (i1 + 1 == i2 || j1 == j2) {
    //align using quadratic space alignment
    var subM = new Array();
    var subN = new Array();

    for (var i = i1 + 1; i <= i2; i++) {
      subM.push(this.M[i - 1]);
    }

    for (var j = j1 + 1; j <= j2; j++) {
      subN.push(this.N[j - 1]);
    }

    var alignment = new AlignPairQuad();

    subScoreSet = new ScoreSet();
    if (j1 == j2) {
      if (j1 == 0) {
        subScoreSet.setScoreSetParam(
          this.scoreSet.scoringMatrix,
          this.scoreSet.beginGap,
          this.scoreSet.beginGap,
          this.scoreSet.beginGap
        );
      } else if (j1 == this.N.length) {
        subScoreSet.setScoreSetParam(
          this.scoreSet.scoringMatrix,
          this.scoreSet.endGap,
          this.scoreSet.endGap,
          this.scoreSet.endGap
        );
      } else {
        subScoreSet.setScoreSetParam(
          this.scoreSet.scoringMatrix,
          this.scoreSet.gap,
          this.scoreSet.gap,
          this.scoreSet.gap
        );
      }
    } else {
      subScoreSet.setScoreSetParam(
        this.scoreSet.scoringMatrix,
        this.scoreSet.gap,
        this.scoreSet.beginGap,
        this.scoreSet.endGap
      );
      subScoreSet.useBeginGapTop = false;
      subScoreSet.useBeginGapLeft = false;
      subScoreSet.useEndGapBottom = false;
      subScoreSet.useEndGapRight = false;

      if (i1 == 0) {
        subScoreSet.useBeginGapTop = true;
      }

      if (j1 == 0) {
        subScoreSet.useBeginGapLeft = true;
      }

      if (j2 == this.N.length) {
        subScoreSet.useEndGapRight = true;
      }

      if (i2 == this.M.length) {
        subScoreSet.useEndGapBottom = true;
      }
    }

    alignment.initializeMatrix(subM, subN, subScoreSet);
    alignment.fillMatrix();
    alignment.align();
    //alignment.dumpMatrix();
    this.alignedM.push(alignment.getAlignedM());
    this.alignedN.push(alignment.getAlignedN());

    this.score = this.score + alignment.score;
  } else {
    var middle = Math.floor((i1 + i2) / 2);

    //linear-space computation of alignment score to middle row
    //forward pass

    //gaps along top

    this.Sn[j1] = 0;

    if (i1 == 0) {
      for (var j = j1 + 1; j <= j2; j++) {
        this.Sn[j] = this.Sn[j - 1] - this.scoreSet.beginGap;
      }
    } else {
      for (var j = j1 + 1; j <= j2; j++) {
        this.Sn[j] = this.Sn[j - 1] - this.scoreSet.gap;
      }
    }

    //now continue down rows to middle row
    var diag;
    var left;
    //for (var i = i1 + 1; i <= i2; i++) {
    for (var i = i1 + 1; i <= middle; i++) {
      diag = this.Sn[j1];
      left;
      if (j1 == 0) {
        left = this.Sn[j1] - this.scoreSet.beginGap;
      } else {
        left = this.Sn[j1] - this.scoreSet.gap;
      }

      this.Sn[j1] = left;

      //we need three values to set the score: diag, left, and above to fill in the row
      for (var j = j1 + 1; j <= j2; j++) {
        //above will be in the this.Sn array, which is holding a mixture of the previous row and the new row
        //var above = this.Sn[j];

        //pick max of three and store in next left
        if (j == this.N.length && i == this.M.length) {
          left = Math.max(
            this.Sn[j] - this.scoreSet.endGap,
            Math.max(
              left - this.scoreSet.endGap,
              diag + this.scoreSet.getScore(this.M[i - 1], this.N[j - 1])
            )
          );
        } else if (i == this.M.length) {
          left = Math.max(
            this.Sn[j] - this.scoreSet.gap,
            Math.max(
              left - this.scoreSet.endGap,
              diag + this.scoreSet.getScore(this.M[i - 1], this.N[j - 1])
            )
          );
        } else if (j == this.N.length) {
          left = Math.max(
            this.Sn[j] - this.scoreSet.endGap,
            Math.max(
              left - this.scoreSet.gap,
              diag + this.scoreSet.getScore(this.M[i - 1], this.N[j - 1])
            )
          );
        } else {
          left = Math.max(
            this.Sn[j] - this.scoreSet.gap,
            Math.max(
              left - this.scoreSet.gap,
              diag + this.scoreSet.getScore(this.M[i - 1], this.N[j - 1])
            )
          );
        }
        diag = this.Sn[j];

        //prepares this.Sn for use in next iteration of i loop
        this.Sn[j] = left;
      }
    }

    //linear-space computation of alignment score to middle row
    //reverse pass

    //gaps along bottom

    this.Sp[j2] = 0;

    if (i2 == this.M.length) {
      for (var j = j2 - 1; j >= j1; j--) {
        this.Sp[j] = this.Sp[j + 1] - this.scoreSet.endGap;
      }
    } else {
      for (var j = j2 - 1; j >= j1; j--) {
        this.Sp[j] = this.Sp[j + 1] - this.scoreSet.gap;
      }
    }

    //now continue up rows to middle row
    var right;
    //for (var i = i2 - 1; i >= i1; i--) {
    for (var i = i2 - 1; i >= middle; i--) {
      diag = this.Sp[j2];
      if (j2 == this.N.length) {
        right = this.Sp[j2] - this.scoreSet.endGap;
      } else {
        right = this.Sp[j2] - this.scoreSet.gap;
      }

      this.Sp[j2] = right;

      //we need three values to set the score: diag, right, and below to fill in the row
      for (var j = j2 - 1; j >= j1; j--) {
        //below will be in the this.Sp array, which is holding a mixture of the previous row and the new row
        //var below = this.Sp[j];

        //pick max of three and store in next right
        if (j == 0 && i == 0) {
          right = Math.max(
            this.Sp[j] - this.scoreSet.beginGap,
            Math.max(
              right - this.scoreSet.beginGap,
              diag +
                this.scoreSet.getScore(this.M[i + 1 - 1], this.N[j + 1 - 1])
            )
          );
        } else if (j == 0) {
          right = Math.max(
            this.Sp[j] - this.scoreSet.beginGap,
            Math.max(
              right - this.scoreSet.gap,
              diag +
                this.scoreSet.getScore(this.M[i + 1 - 1], this.N[j + 1 - 1])
            )
          );
        } else if (i == 0) {
          right = Math.max(
            this.Sp[j] - this.scoreSet.gap,
            Math.max(
              right - this.scoreSet.beginGap,
              diag +
                this.scoreSet.getScore(this.M[i + 1 - 1], this.N[j + 1 - 1])
            )
          );
        } else {
          right = Math.max(
            this.Sp[j] - this.scoreSet.gap,
            Math.max(
              right - this.scoreSet.gap,
              diag +
                this.scoreSet.getScore(this.M[i + 1 - 1], this.N[j + 1 - 1])
            )
          );
        }
        diag = this.Sp[j];
        this.Sp[j] = right;
      }
    }

    //now find the value of j that maximizes this.Sn[j] + this.Sp[j]
    //this point will be in the maximum scoring path in the final alignment.
    //once we have this point we can divide the problem into two new problems,

    var maxValue = this.Sn[j1] + this.Sp[j1];
    var maxJ = j1;

    for (var j = j1 + 1; j <= j2; j++) {
      if (this.Sn[j] + this.Sp[j] >= maxValue) {
        maxValue = this.Sn[j] + this.Sp[j];
        maxJ = j;
      }
    }

    this.path(i1, j1, middle, maxJ);
    this.path(middle, maxJ, i2, j2);
  }
}

//AlignPairLinear class getAlignedM() method
function getAlignedM() {
  return this.alignedM.join("");
}

//AlignPairLinear class getAlignedN() method
function getAlignedN() {
  return this.alignedN.join("");
}

//AlignPairLinear class setAlignParam method
function setAlignParam(M, N, scoreSet) {
  this.M = M;
  this.N = N;
  this.alignedM = new Array();
  this.alignedN = new Array();
  this.scoreSet = scoreSet;
  this.Sn = new Array(this.N.length);
  this.Sp = new Array(this.N.length);
  this.score = 0;
}

//AlignPairLinear class
function AlignPairLinear() {
  this.M;
  this.N;
  this.alignedM;
  this.alignedN;
  this.scoreSet;
  this.Sn;
  this.Sp;
  this.score;
}

//create and throw away a prototype object
new AlignPairLinear();

//define object methods
AlignPairLinear.prototype.align = align;
AlignPairLinear.prototype.path = path;
AlignPairLinear.prototype.setAlignParam = setAlignParam;
AlignPairLinear.prototype.getAlignedM = getAlignedM;
AlignPairLinear.prototype.getAlignedN = getAlignedN;
