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

//This class should be used for small alignments,
//since it uses O(nm) memory, where n and m are the sequence lengths.
//For larger alignments use the linear space algorithm implemented
//in align_pair_linear.js

//To use this class: (see pairwise_dna.js for example)
//var alignment = new AlignPairQuad();
//alignment.initializeMatrix(sequenceArrayM, sequenceArrayN, scoreSet);
//alignment.fillMatrix();
//alignment.align();
//var alignedSequenceStringM = alignment.getAlignedM();
//var alignedSequenceStringN = alignment.getAlignedN();

//------------------------------------ Node class
//Node class
function Node() {
  this.value;
  this.tracebackI;
  this.tracebackJ;
}
//------------------------------------

//------------------------------------ AlignPairQuad class
//AlignPairQuad class initializeMatrix method
function initializeMatrix(sequenceOne, sequenceTwo, scoreSet) {
  this.scoreSet = scoreSet;

  this.M = sequenceOne;
  this.N = sequenceTwo;
  this.score = 0;

  //create an two-dimensional array of nodes
  this.nodes = new Array(this.M.length + 1);

  //row i
  for (var i = 0; i < this.nodes.length; i++) {
    this.nodes[i] = new Array(this.N.length + 1);
    //column j
    for (var j = 0; j < this.nodes[i].length; j++) {
      this.nodes[i][j] = new Node();
    }
  }

  this.nodes[0][0].value = 0;

  //i rows
  for (var i = 1; i < this.nodes.length; i++) {
    if (this.scoreSet.useBeginGapLeft) {
      this.nodes[i][0].value =
        this.nodes[i - 1][0].value - this.scoreSet.beginGap;
    } else {
      this.nodes[i][0].value = this.nodes[i - 1][0].value - this.scoreSet.gap;
    }
    this.nodes[i][0].tracebackI = i - 1;
    this.nodes[i][0].tracebackJ = 0;
  }

  //j columns
  for (var j = 1; j < this.nodes[0].length; j++) {
    if (this.scoreSet.useBeginGapTop) {
      this.nodes[0][j].value =
        this.nodes[0][j - 1].value - this.scoreSet.beginGap;
    } else {
      this.nodes[0][j].value = this.nodes[0][j - 1].value - this.scoreSet.gap;
    }
    this.nodes[0][j].tracebackI = 0;
    this.nodes[0][j].tracebackJ = j - 1;
  }
}

//AlignPairQuad class dumpMatrix method
function dumpMatrix() {
  outputWindow.document.write(
    "Dynamic programming matrix i=" +
      this.nodes.length +
      " and j=" +
      this.nodes[0].length
  );
  outputWindow.document.write("\n");
  for (var i = 0; i < this.nodes.length; i++) {
    for (var j = 0; j < this.nodes[i].length; j++) {
      var traceI = this.nodes[i][j].tracebackI;
      var traceJ = this.nodes[i][j].tracebackJ;

      if (traceI == undefined) {
        traceI = "u";
      }
      if (traceJ == undefined) {
        traceJ = "u";
      }
      var output =
        "(" +
        i +
        "," +
        j +
        ")[" +
        traceI +
        "," +
        traceJ +
        "]=" +
        this.nodes[i][j].value;
      outputWindow.document.write(rightNum(output, "", 20, " "));
    }
    outputWindow.document.write("\n");
  }
  outputWindow.document.write("\n");
}

//AlignPairQuad class fillMatrix method
function fillMatrix() {
  //i rows
  for (var i = 1; i < this.nodes.length; i++) {
    //j columns
    for (var j = 1; j < this.nodes[0].length; j++) {
      var a;
      var b;
      var c;

      //handle end gaps here

      if (i == this.nodes.length - 1 && j == this.nodes[0].length - 1) {
        if (this.scoreSet.useEndGapRight) {
          a = this.nodes[i - 1][j].value - this.scoreSet.endGap;
        } else {
          a = this.nodes[i - 1][j].value - this.scoreSet.gap;
        }

        if (this.scoreSet.useEndGapBottom) {
          b = this.nodes[i][j - 1].value - this.scoreSet.endGap;
        } else {
          b = this.nodes[i][j - 1].value - this.scoreSet.gap;
        }
      } else if (i == this.nodes.length - 1) {
        a = this.nodes[i - 1][j].value - this.scoreSet.gap;
        if (this.scoreSet.useEndGapBottom) {
          b = this.nodes[i][j - 1].value - this.scoreSet.endGap;
        } else {
          b = this.nodes[i][j - 1].value - this.scoreSet.gap;
        }
      } else if (j == this.nodes[0].length - 1) {
        if (this.scoreSet.useEndGapRight) {
          a = this.nodes[i - 1][j].value - this.scoreSet.endGap;
        } else {
          a = this.nodes[i - 1][j].value - this.scoreSet.gap;
        }
        b = this.nodes[i][j - 1].value - this.scoreSet.gap;
      } else {
        a = this.nodes[i - 1][j].value - this.scoreSet.gap;
        b = this.nodes[i][j - 1].value - this.scoreSet.gap;
      }

      c =
        this.nodes[i - 1][j - 1].value +
        this.scoreSet.getScore(this.M[i - 1], this.N[j - 1]);

      if (a >= b && a >= c) {
        this.nodes[i][j].value = a;
        this.nodes[i][j].tracebackI = i - 1;
        this.nodes[i][j].tracebackJ = j;
      } else if (b >= c && b >= a) {
        this.nodes[i][j].value = b;
        this.nodes[i][j].tracebackI = i;
        this.nodes[i][j].tracebackJ = j - 1;
      } else {
        this.nodes[i][j].value = c;
        this.nodes[i][j].tracebackI = i - 1;
        this.nodes[i][j].tracebackJ = j - 1;
      }
    }
  }
  this.score = this.nodes[this.nodes.length - 1][
    this.nodes[0].length - 1
  ].value;
}

//AlignPairQuad class align() method
function align() {
  this.alignedM = new Array();
  this.alignedN = new Array();

  var currentI = this.nodes.length - 1;
  var currentJ = this.nodes[0].length - 1;

  var currentNode = this.nodes[this.nodes.length - 1][this.nodes[0].length - 1];

  while (
    currentNode.tracebackI != undefined &&
    currentNode.tracebackJ != undefined
  ) {
    if (
      currentNode.tracebackI == currentI - 1 &&
      currentNode.tracebackJ == currentJ - 1
    ) {
      this.alignedM.push(this.M.pop());
      this.alignedN.push(this.N.pop());
    }
    // edited here .-- instead of - because of codon width
    else if (currentNode.tracebackJ == currentJ - 1) {
      this.alignedM.push(".--");
      this.alignedN.push(this.N.pop());
    } else {
      this.alignedM.push(this.M.pop());
      this.alignedN.push(".--");
    }

    currentI = currentNode.tracebackI;
    currentJ = currentNode.tracebackJ;

    currentNode = this.nodes[currentNode.tracebackI][currentNode.tracebackJ];
  }

  this.alignedM = this.alignedM.reverse();
  this.alignedN = this.alignedN.reverse();
}

//AlignPairQuad class getAlignedM() method
function getAlignedM() {
  return this.alignedM.join("");
}

//AlignPairQuad class getAlignedN() method
function getAlignedN() {
  return this.alignedN.join("");
}

//AlignPairQuad class
function AlignPairQuad() {
  this.M;
  this.N;
  this.scoreSet;
  this.nodes;
  this.alignedM;
  this.alignedN;
  this.score;
}

//create and throw away a prototype object
new AlignPairQuad();

//define object methods
AlignPairQuad.prototype.initializeMatrix = initializeMatrix;
AlignPairQuad.prototype.fillMatrix = fillMatrix;
AlignPairQuad.prototype.align = align;
AlignPairQuad.prototype.getAlignedM = getAlignedM;
AlignPairQuad.prototype.getAlignedN = getAlignedN;
AlignPairQuad.prototype.dumpMatrix = dumpMatrix;
