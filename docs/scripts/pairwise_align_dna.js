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

function pairwiseAlignDna(theDocument) {
  //var MATCH_SCORE = 2;
  //var MISMATCH_SCORE = -1;
  //var GAP_PENALTY = 2;

  //var BEGIN_GAP_PENALTY = 0;
  //var END_GAP_PENALTY = 0;

  var newDnaOne = "";
  var titleOne = "";

  var newDnaTwo = "";
  var titleTwo = "";

  var maxInput = 20000;

  if (testScript() == false) {
    return false;
  }

  if (
    checkFormElement(theDocument.forms[0].elements[0]) == false ||
    checkSequenceLength(theDocument.forms[0].elements[0].value, maxInput) ==
      false ||
    checkFormElement(theDocument.forms[0].elements[1]) == false ||
    checkSequenceLength(theDocument.forms[0].elements[1].value, maxInput) ==
      false
  ) {
    return false;
  }

  var MATCH_SCORE = parseInt(
    theDocument.forms[0].elements[5].options[
      theDocument.forms[0].elements[5].selectedIndex
    ].value
  );
  var MISMATCH_SCORE = parseInt(
    theDocument.forms[0].elements[6].options[
      theDocument.forms[0].elements[6].selectedIndex
    ].value
  );
  var BEGIN_GAP_PENALTY = parseInt(
    theDocument.forms[0].elements[7].options[
      theDocument.forms[0].elements[7].selectedIndex
    ].value
  );
  var GAP_PENALTY = parseInt(
    theDocument.forms[0].elements[8].options[
      theDocument.forms[0].elements[8].selectedIndex
    ].value
  );
  var END_GAP_PENALTY = parseInt(
    theDocument.forms[0].elements[9].options[
      theDocument.forms[0].elements[9].selectedIndex
    ].value
  );

  openWindow("Pairwise Align DNA");
  openPre();

  newDnaOne = getSequenceFromFasta(theDocument.forms[0].elements[0].value);
  newDnaOne = removeNonDna(newDnaOne);
  titleOne = getTitleFromFasta(theDocument.forms[0].elements[0].value);

  newDnaTwo = getSequenceFromFasta(theDocument.forms[0].elements[1].value);
  newDnaTwo = removeNonDna(newDnaTwo);
  titleTwo = getTitleFromFasta(theDocument.forms[0].elements[1].value);

  outputWindow.document.write(
    getPairwiseAlignTitle(titleOne, newDnaOne, titleTwo, newDnaTwo)
  );

  //change to arrays for pass by reference, so that large sequence isn't copied
  if (newDnaOne.search(/./) != -1) {
    newDnaOne = newDnaOne.match(/./g);
  }

  if (newDnaTwo.search(/./) != -1) {
    newDnaTwo = newDnaTwo.match(/./g);
  }

  pairwiseDna(
    titleOne,
    newDnaOne,
    titleTwo,
    newDnaTwo,
    MATCH_SCORE,
    MISMATCH_SCORE,
    GAP_PENALTY,
    BEGIN_GAP_PENALTY,
    END_GAP_PENALTY
  );
  closePre();
  closeWindow();
  return true;
}

function pairwiseDna(
  titleOne,
  newDnaOne,
  titleTwo,
  newDnaTwo,
  matchScore,
  mismatchScore,
  gapPenalty,
  beginGapPenalty,
  endGapPenalty
) {
  //can use one or both.
  //can compare scores (should be identical)
  var useLinearSpace = true;
  var useQuadraticSpace = false;

  var matrix = new Identity();
  matrix.setMatch(matchScore);
  matrix.setMismatch(mismatchScore);

  var scoreSet = new ScoreSet();
  scoreSet.setScoreSetParam(matrix, gapPenalty, beginGapPenalty, endGapPenalty);

  var alignment;

  if (useLinearSpace) {
    alignment = new AlignPairLinear();
    alignment.setAlignParam(newDnaOne, newDnaTwo, scoreSet);
    alignment.align();

    outputWindow.document.write(">" + titleOne + "\n");
    outputWindow.document.write(addReturns(alignment.getAlignedM()));
    outputWindow.document.write("\n");
    outputWindow.document.write("\n");
    outputWindow.document.write(">" + titleTwo + "\n");
    outputWindow.document.write(addReturns(alignment.getAlignedN()));
    outputWindow.document.write("\n\n");
    outputWindow.document.write("Alignment score: " + alignment.score + "\n\n");
  }

  if (useQuadraticSpace) {
    alignment = new AlignPairQuad();
    alignment.initializeMatrix(newDnaOne, newDnaTwo, scoreSet);
    alignment.fillMatrix();
    //alignment.dumpMatrix();
    alignment.align();

    outputWindow.document.write(">" + titleOne + "\n");
    outputWindow.document.write(addReturns(alignment.getAlignedM()));
    outputWindow.document.write("\n");
    outputWindow.document.write("\n");
    outputWindow.document.write(">" + titleTwo + "\n");
    outputWindow.document.write(addReturns(alignment.getAlignedN()));
    outputWindow.document.write("\n\n");
    outputWindow.document.write("Alignment score: " + alignment.score + "\n\n");
  }
}

//------------------------------------ ScoreSet class

//ScoreSet getScore
function getScore(r1, r2) {
  return this.scoringMatrix.scoringMatrix_getScore(r1, r2);
}

//ScoreSet setScoreSetParam
function setScoreSetParam(
  scoringMatrix,
  gapPenalty,
  beginGapPenalty,
  endGapPenalty
) {
  this.scoringMatrix = scoringMatrix;
  this.gap = gapPenalty;
  this.beginGap = beginGapPenalty;
  this.endGap = endGapPenalty;
}

//ScoreSet class
function ScoreSet() {
  this.scoringMatrix;
  this.gap;
  this.beginGap;
  this.endGap;
  this.useBeginGapTop = true;
  this.useBeginGapLeft = true;
  this.useEndGapBottom = true;
  this.useEndGapRight = true;
}

//create and throw away a prototype object
new ScoreSet();

//define object methods
ScoreSet.prototype.getScore = getScore;
ScoreSet.prototype.setScoreSetParam = setScoreSetParam;

//------------------------------------

//------------------------------------ ScoringMatrix abstract class
//ScoringMatrix getScore method
function scoringMatrix_getScore(r1, r2) {
  r1 = r1.toLowerCase();
  r2 = r2.toLowerCase();
  if (r1 == r2) {
    return this.match;
  } else {
    return this.mismatch;
  }
}

//ScoringMatrix class
function ScoringMatrix() {
  this.mismatch;
  this.match;
}

//create and throw away a prototype object
new ScoringMatrix();

//define object methods
ScoringMatrix.prototype.scoringMatrix_getScore = scoringMatrix_getScore;

//------------------------------------ Identity class extends ScoringMatrix Class
//Identity class setMismatch method
function setMismatch(mismatchScore) {
  this.mismatch = mismatchScore;
}

//Identity class setMatch method
function setMatch(matchScore) {
  this.match = matchScore;
}

//Identity class
function Identity() {}

Identity.prototype = new ScoringMatrix();
Identity.prototype.setMismatch = setMismatch;
Identity.prototype.setMatch = setMatch;
