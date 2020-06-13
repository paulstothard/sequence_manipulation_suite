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

function randomCodingDna(theDocument) {
  var maxInput = 4000000;

  if (testScript() == false) {
    return false;
  }

  var enteredNumber = theDocument.forms[0].elements[0].value.replace(
    /[^\d]/g,
    ""
  );
  if (
    checkFormElement(theDocument.forms[0].elements[0]) == false ||
    verifyMaxDigits(enteredNumber, maxInput) == false
  ) {
    return false;
  }

  var seqNum = parseInt(
    theDocument.forms[0].elements[5].options[
      theDocument.forms[0].elements[5].selectedIndex
    ].value
  );
  var geneticCode = getGeneticCodeString(
    theDocument.forms[0].elements[4].options[
      theDocument.forms[0].elements[4].selectedIndex
    ].value
  );
  geneticCode = geneticCode.split(/,/);
  var codonListSet = new GeneticCode(geneticCode);
  codonListSet.parseGeneticCodeArray();

  openWindow("Random Coding DNA");
  openPre();
  for (var i = 1; i <= seqNum; i++) {
    outputWindow.document.write(
      "&gt;" +
        "random coding sequence " +
        i +
        " consisting of " +
        enteredNumber * 3 +
        " bases.\n"
    );
    writeRandomCodingDna(
      codonListSet.startCodons,
      codonListSet.stopCodons,
      codonListSet.codingCodons,
      enteredNumber
    );
    outputWindow.document.write("\n");
  }
  closePre();
  closeWindow();
  return true;
}

function writeRandomCodingDna(
  startCodons,
  stopCodons,
  codingCodons,
  lengthInCodons
) {
  var sequence = "";
  var tempNum = 0;
  var tempChar = "";

  lengthInCodons = parseInt(lengthInCodons);

  //add start codon
  if (lengthInCodons > 0) {
    tempNum = Math.round(Math.random() * startCodons.length);
    if (tempNum == startCodons.length) {
      tempNum = 0;
    }
    sequence = sequence + startCodons[tempNum];
  }

  //add coding codons
  for (var j = 0; j < lengthInCodons - 2; j++) {
    tempNum = Math.round(Math.random() * codingCodons.length);
    if (tempNum == codingCodons.length) {
      tempNum = 0;
    }
    sequence = sequence + codingCodons[tempNum];
    if (sequence.length == 60) {
      outputWindow.document.write(sequence + "\n");
      sequence = "";
    }
  }

  //add stop codon
  if (lengthInCodons > 1) {
    tempNum = Math.round(Math.random() * stopCodons.length);
    if (tempNum == stopCodons.length) {
      tempNum = 0;
    }
    sequence = sequence + stopCodons[tempNum];
  }

  outputWindow.document.write(sequence + "\n");
  return true;
}

//class GeneticCode method parseGeneticCodeArray()
function parseGeneticCodeArray() {
  var codonSequence =
    "gggggaggtggcgaggaagatgacgtggtagttgtcgcggcagctgccaggagaagtagcaagaaaaataacatgataattatcacgacaactacctggtgatgttgctagtaatattacttgttatttttctcgtcatcttcccggcgacgtcgccagcaacatcacctgctacttctcccgccacctccc";
  var proteinSequence;
  var geneticCodeMatchExp = getGeneticCodeMatchExp(this.geneticCodeArray);
  var geneticCodeMatchResult = getGeneticCodeMatchResult(this.geneticCodeArray);

  codonSequence = codonSequence.replace(/(...)/g, function (
    str,
    p1,
    offset,
    s
  ) {
    return " " + p1 + " ";
  });

  var codonSequenceCopy = codonSequence;

  for (var i = 0; i < geneticCodeMatchExp.length; i++) {
    codonSequence = codonSequence.replace(
      geneticCodeMatchExp[i],
      geneticCodeMatchResult[i]
    );
  }
  var codonArray = codonSequenceCopy.split(/\s+/);

  codonSequence = codonSequence.replace(/\*/g, "Z");
  var proteinArray = codonSequence.split(/\s+/);

  for (var i = 0; i < codonArray.length; i++) {
    //on some systems there will be empty items because of split function
    if (proteinArray[i] == "" && codonArray[i] == "") {
      continue;
    }
    if (proteinArray[i].toLowerCase() == "z") {
      this.stopCodons.push(codonArray[i]);
    } else if (proteinArray[i].toLowerCase() == "m") {
      this.startCodons.push(codonArray[i]);
      this.codingCodons.push(codonArray[i]);
    } else {
      this.codingCodons.push(codonArray[i]);
    }
  }
}

//class GeneticCode
function GeneticCode(geneticCodeArray) {
  this.geneticCodeArray = geneticCodeArray;
  this.startCodons = new Array();
  this.stopCodons = new Array();
  //coding will include starts
  this.codingCodons = new Array();
}

//create and throw away a prototype object
new GeneticCode();

// define object methods
GeneticCode.prototype.parseGeneticCodeArray = parseGeneticCodeArray;
