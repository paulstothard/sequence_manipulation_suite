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

//This class should be used when searching for short motifs,
//since it uses O(nm) memory, where n and m are the sequence lengths.

//To use this class: (see fuzzy_search_dna.js for example)
//var fuzzySearch = new FuzzySearch();
//fuzzySearch.initializeMatrix(sequenceArrayM, sequenceArrayN, scoreSet);
//fuzzySearch.search();
//var arrayOfHits = fuzzySearch.getHits();

//------------------------------------ Node class
//Node class
function Node() {
  this.value;
  this.tracebackI;
  this.tracebackJ;
  this.alreadyMatched = false;
}
//------------------------------------

//------------------------------------ FuzzySearch class
//FuzzySearch class initializeMatrix method
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
    this.nodes[i][0].value = this.nodes[i - 1][0].value;
    this.nodes[i][0].tracebackI = i - 1;
    this.nodes[i][0].tracebackJ = 0;
  }

  //j columns
  for (var j = 1; j < this.nodes[0].length; j++) {
    this.nodes[0][j].value = this.nodes[0][j - 1].value;
    this.nodes[0][j].tracebackI = 0;
    this.nodes[0][j].tracebackJ = j - 1;
  }
}

//FuzzySearch class dumpMatrix method
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

//FuzzySearch class updateMatrix method
function updateMatrix() {
  //i rows
  for (var i = 1; i < this.nodes.length; i++) {
    //j columns
    for (var j = 1; j < this.nodes[0].length; j++) {
      var a;
      var b;
      var c;

      if (this.nodes[i][j].alreadyMatched) {
        a = 0;
        b = 0;
        c = 0;
      } else if (i == this.nodes.length - 1 && j == this.nodes[0].length - 1) {
        a = this.nodes[i - 1][j].value;
        b = this.nodes[i][j - 1].value;
        c =
          this.nodes[i - 1][j - 1].value +
          this.scoreSet.getScore(this.M[i - 1], this.N[j - 1]);
      } else if (i == this.nodes.length - 1) {
        a = this.nodes[i - 1][j].value - this.scoreSet.gap;
        b = this.nodes[i][j - 1].value;
        c =
          this.nodes[i - 1][j - 1].value +
          this.scoreSet.getScore(this.M[i - 1], this.N[j - 1]);
      } else if (j == this.nodes[0].length - 1) {
        a = this.nodes[i - 1][j].value;
        b = this.nodes[i][j - 1].value - this.scoreSet.gap;
        c =
          this.nodes[i - 1][j - 1].value +
          this.scoreSet.getScore(this.M[i - 1], this.N[j - 1]);
      } else {
        a = this.nodes[i - 1][j].value - this.scoreSet.gap;
        b = this.nodes[i][j - 1].value - this.scoreSet.gap;
        c =
          this.nodes[i - 1][j - 1].value +
          this.scoreSet.getScore(this.M[i - 1], this.N[j - 1]);
      }

      if (a > b && a > c) {
        this.nodes[i][j].value = a;
        this.nodes[i][j].tracebackI = i - 1;
        this.nodes[i][j].tracebackJ = j;
      } else if (b > c && b > a) {
        this.nodes[i][j].value = b;
        this.nodes[i][j].tracebackI = i;
        this.nodes[i][j].tracebackJ = j - 1;
      } else {
        this.nodes[i][j].value = c;
        this.nodes[i][j].tracebackI = i - 1;
        this.nodes[i][j].tracebackJ = j - 1;
      }
      if (this.nodes[i][j].value < 0) {
        this.nodes[i][j].value = 0;
        this.nodes[i][j].tracebackI = undefined;
        this.nodes[i][j].tracebackJ = undefined;
      }
    }
  }
  this.score = this.nodes[this.nodes.length - 1][
    this.nodes[0].length - 1
  ].value;
}

//FuzzySearch class search() method
function search() {
  this.hits = new Array();
  var hitCount = 0;

  while (hitCount < this.scoreSet.hits) {
    this.updateMatrix();
    //this.dumpMatrix();

    //start at max node on
    var maxNodeValue = 0;
    var maxNodeI = undefined;
    var maxNodeJ = undefined;

    for (var i = 1; i < this.nodes.length; i++) {
      for (var j = 1; j < this.nodes[0].length; j++) {
        if (this.nodes[i][j].value > maxNodeValue) {
          maxNodeValue = this.nodes[i][j].value;
          maxNodeI = i;
          maxNodeJ = j;
        }
      }
    }

    if (maxNodeValue == 0) {
      break;
    }

    var currentI = maxNodeI;
    var currentJ = maxNodeJ;
    var currentNode = this.nodes[maxNodeI][maxNodeJ];

    //do trace back

    var alignedM = new Array();
    var alignedN = new Array();
    var score = currentNode.value;
    var endM = undefined;
    var endN = undefined;

    while (
      currentNode.tracebackI != undefined &&
      currentNode.tracebackJ != undefined
    ) {
      //don't show begin gaps
      if (currentI == 0 || currentJ == 0) {
        break;
      }

      if (
        currentNode.tracebackI == currentI - 1 &&
        currentNode.tracebackJ == currentJ - 1
      ) {
        if (endM == undefined) {
          endM = currentI;
          endN = currentJ;
        }
        alignedM.push(this.M[currentI - 1]);
        alignedN.push(this.N[currentJ - 1]);
      } else if (currentNode.tracebackJ == currentJ - 1) {
        //don't show end gaps
        if (endM != undefined) {
          alignedM.push("-");
          alignedN.push(this.N[currentJ - 1]);
        }
      } else {
        //don't show end gaps
        if (endM != undefined) {
          alignedM.push(this.M[currentI - 1]);
          alignedN.push("-");
        }
      }

      this.nodes[currentI][currentJ].value = 0;
      this.nodes[currentI][currentJ].alreadyMatched = true;

      currentI = currentNode.tracebackI;
      currentJ = currentNode.tracebackJ;

      currentNode = this.nodes[currentNode.tracebackI][currentNode.tracebackJ];
    }

    alignedM = alignedM.reverse();
    alignedN = alignedN.reverse();

    this.hits.push(
      new Hit(
        alignedM.join(""),
        alignedN.join(""),
        score,
        currentI + 1,
        endM,
        currentJ + 1,
        endN
      )
    );

    hitCount++;
  }
}

//FuzzySearch class getHits() method
function getHits() {
  return this.hits;
}

//FuzzySearch class
function FuzzySearch() {
  this.M;
  this.N;
  this.scoreSet;
  this.nodes;
  this.hits;
}

//create and throw away a prototype object
new FuzzySearch();

//define object methods
FuzzySearch.prototype.initializeMatrix = initializeMatrix;
FuzzySearch.prototype.updateMatrix = updateMatrix;
FuzzySearch.prototype.search = search;
FuzzySearch.prototype.getHits = getHits;
FuzzySearch.prototype.dumpMatrix = dumpMatrix;
//------------------------------------

//------------------------------------ Hit class
//Hit class
function Hit(sequenceM, sequenceN, score, startM, endM, startN, endN) {
  this.sequenceM = sequenceM;
  this.sequenceN = sequenceN;
  this.score = score;
  this.startM = startM;
  this.endM = endM;
  this.startN = startN;
  this.endN = endN;
}
//------------------------------------
