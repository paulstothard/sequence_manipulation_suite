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

function fuzzySearchDna(theDocument) {
  //var MATCH_SCORE = 2;
  //var MISMATCH_SCORE = -1;
  //var GAP_PENALTY = 2;
  //var HITS = 8;

  var targetSequence = "";
  var targetTitle = "";

  var querySequence = "";
  var queryTitle = "";

  var maxTarget = 2000000;
  var maxQuery = 30;

  if (testScript() == false) {
    return false;
  }

  if (
    checkFormElement(theDocument.forms[0].elements[0]) == false ||
    checkSequenceLength(theDocument.forms[0].elements[0].value, maxTarget) ==
      false ||
    checkFormElement(theDocument.forms[0].elements[1]) == false ||
    checkSequenceLength(theDocument.forms[0].elements[1].value, maxQuery) ==
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
  var GAP_PENALTY = parseInt(
    theDocument.forms[0].elements[7].options[
      theDocument.forms[0].elements[7].selectedIndex
    ].value
  );
  var HITS = parseInt(
    theDocument.forms[0].elements[8].options[
      theDocument.forms[0].elements[8].selectedIndex
    ].value
  );

  openWindow("Fuzzy Search DNA");
  openPre();

  targetSequence = getSequenceFromFasta(theDocument.forms[0].elements[0].value);
  targetSequence = removeNonDna(targetSequence);
  targetTitle = getTitleFromFasta(theDocument.forms[0].elements[0].value);

  querySequence = getSequenceFromFasta(theDocument.forms[0].elements[1].value);
  querySequence = removeNonDna(querySequence);
  queryTitle = "query";

  outputWindow.document.write(
    getFuzzySearchTitle(targetTitle, targetSequence, queryTitle, querySequence)
  );

  //change to arrays for pass by reference, so that large sequence isn't copied
  if (targetSequence.search(/./) != -1) {
    targetSequence = targetSequence.match(/./g);
  }

  if (querySequence.search(/./) != -1) {
    querySequence = querySequence.match(/./g);
  }

  if (targetSequence.length == 0) {
    alert("The sequence contains no DNA bases.");
    return false;
  }

  if (querySequence.length == 0) {
    alert("The query sequence contains no DNA bases.");
    return false;
  }

  _fuzzySearchDna(
    queryTitle,
    querySequence,
    targetTitle,
    targetSequence,
    MATCH_SCORE,
    MISMATCH_SCORE,
    GAP_PENALTY,
    HITS
  );
  closePre();
  closeWindow();
  return true;
}

function _fuzzySearchDna(
  queryTitle,
  querySequence,
  targetTitle,
  targetSequence,
  matchScore,
  mismatchScore,
  gapPenalty,
  hits
) {
  var matrix = new Identity();
  matrix.setMatch(matchScore);
  matrix.setMismatch(mismatchScore);

  var scoreSet = new ScoreSet();
  scoreSet.setScoreSetParam(matrix, gapPenalty, hits);

  var fuzzySearch = new FuzzySearch();
  fuzzySearch.initializeMatrix(querySequence, targetSequence, scoreSet);
  fuzzySearch.search();
  var hits = fuzzySearch.getHits();

  if (hits.length > 0) {
    for (var i = 0; i < hits.length; i++) {
      outputWindow.document.write(
        ">" +
          queryTitle +
          " from " +
          hits[i].startM +
          " to " +
          hits[i].endM +
          "\n"
      );
      outputWindow.document.write(hits[i].sequenceM + "\n");
      outputWindow.document.write(
        ">" +
          targetTitle +
          " from " +
          hits[i].startN +
          " to " +
          hits[i].endN +
          "\n"
      );
      outputWindow.document.write(hits[i].sequenceN + "\n");
      outputWindow.document.write("Score: " + hits[i].score + "\n\n");
    }
  } else {
    outputWindow.document.write("No hits were obtained.\n\n");
  }
}

//------------------------------------ ScoreSet class

//ScoreSet getScore
function getScore(r1, r2) {
  return this.scoringMatrix.scoringMatrix_getScore(r1, r2);
}

//ScoreSet setScoreSetParam
function setScoreSetParam(scoringMatrix, gapPenalty, hits) {
  this.scoringMatrix = scoringMatrix;
  this.gap = gapPenalty;
  this.hits = hits;
}

//ScoreSet class
function ScoreSet() {
  this.scoringMatrix;
  this.gap;
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
