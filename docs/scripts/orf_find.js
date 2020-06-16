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

function orfFind(theDocument) {
  var newDna = "";
  var title = "";
  var maxInput = 100000000;

  if (testScript() == false) {
    return false;
  }

  var geneticCode = getGeneticCodeString(
    theDocument.forms[0].elements[8].options[
      theDocument.forms[0].elements[8].selectedIndex
    ].value
  );

  geneticCode = geneticCode.split(/,/);

  var enteredNumber = theDocument.forms[0].elements[7].value.replace(
    /[^\d]/g,
    ""
  );
  if (
    checkFormElement(theDocument.forms[0].elements[0]) == false ||
    checkFormElement(theDocument.forms[0].elements[7]) == false ||
    verifyDigits(enteredNumber) == false ||
    checkSequenceLength(theDocument.forms[0].elements[0].value, maxInput) ==
      false
  ) {
    return false;
  }

  if (checkGeneticCode(geneticCode) == false) {
    return false;
  }

  newDna = getSequenceFromFasta(theDocument.forms[0].elements[0].value);
  title = getTitleFromFasta(theDocument.forms[0].elements[0].value);
  verifyDna(newDna);
  newDna = removeNonDna(newDna);
  openWindow("ORF Finder");
  outputWindow.document.write(getInfoFromTitleAndSequence(title, newDna));
  openPre();

  if (
    theDocument.forms[0].elements[5].options[
      theDocument.forms[0].elements[5].selectedIndex
    ].value == "all"
  ) {
    //rf 1
    writeOrfs(
      newDna,
      geneticCode,
      theDocument.forms[0].elements[4].options[
        theDocument.forms[0].elements[4].selectedIndex
      ].value,
      0,
      theDocument.forms[0].elements[6].options[
        theDocument.forms[0].elements[6].selectedIndex
      ].value,
      enteredNumber
    );
    //rf 2
    writeOrfs(
      newDna,
      geneticCode,
      theDocument.forms[0].elements[4].options[
        theDocument.forms[0].elements[4].selectedIndex
      ].value,
      1,
      theDocument.forms[0].elements[6].options[
        theDocument.forms[0].elements[6].selectedIndex
      ].value,
      enteredNumber
    );
    //rf 3
    writeOrfs(
      newDna,
      geneticCode,
      theDocument.forms[0].elements[4].options[
        theDocument.forms[0].elements[4].selectedIndex
      ].value,
      2,
      theDocument.forms[0].elements[6].options[
        theDocument.forms[0].elements[6].selectedIndex
      ].value,
      enteredNumber
    );
  } else {
    writeOrfs(
      newDna,
      geneticCode,
      theDocument.forms[0].elements[4].options[
        theDocument.forms[0].elements[4].selectedIndex
      ].value,
      theDocument.forms[0].elements[5].options[
        theDocument.forms[0].elements[5].selectedIndex
      ].value,
      theDocument.forms[0].elements[6].options[
        theDocument.forms[0].elements[6].selectedIndex
      ].value,
      enteredNumber
    );
  }

  closePre();
  closeWindow();
  return true;
}

function writeOrfs(
  dnaSequence,
  geneticCode,
  startCodons,
  startPos,
  strand,
  theLength
) {
  var i = 0;
  var k = 0;
  var codon = "";
  var foundStart = false;
  var geneticCodeMatchExp = getGeneticCodeMatchExp(geneticCode);
  var geneticCodeMatchResult = getGeneticCodeMatchResult(geneticCode);
  var proteinLength = 0;
  var foundStop = false;

  var geneticCodeMatchExpStop;
  for (var j = 0; j < geneticCodeMatchExp.length; j++) {
    if (geneticCodeMatchResult[j] == "*") {
      geneticCodeMatchExpStop = geneticCodeMatchExp[j];
      break;
    }
  }

  var startRe = new RegExp(startCodons, "i");
  var sequenceToTranslate;

  startPos = parseInt(startPos);
  var rf = startPos + 1;
  theLength = parseInt(theLength);

  if (strand == "reverse") {
    dnaSequence = reverse(complement(dnaSequence));
  }

  while (i <= dnaSequence.length - 3) {
    for (var i = startPos; i <= dnaSequence.length - 3; i = i + 3) {
      codon = dnaSequence.substring(i, i + 3);
      if (
        startCodons != "any" &&
        foundStart == false &&
        codon.search(startRe) == -1
      ) {
        break;
      }
      foundStart = true;

      if (codon.search(geneticCodeMatchExpStop) != -1) {
        foundStop = true;
      }

      proteinLength++;

      if (foundStop && proteinLength < theLength) {
        break;
      }
      if (
        (foundStop && proteinLength >= theLength) ||
        (i >= dnaSequence.length - 5 && proteinLength >= theLength)
      ) {
        sequenceToTranslate = dnaSequence.substring(startPos, i + 3);

        outputWindow.document.write(
          "&gt;ORF number " +
            (k + 1) +
            " in reading frame " +
            rf +
            " on the " +
            strand +
            " strand extends from base " +
            (startPos + 1) +
            " to base " +
            (i + 3) +
            ".\n"
        );

        outputWindow.document.write(addReturns(sequenceToTranslate) + "\n\n");

        outputWindow.document.write(
          "&gt;Translation of ORF number " +
            (k + 1) +
            " in reading frame " +
            rf +
            " on the " +
            strand +
            " strand.\n"
        );

        sequenceToTranslate = sequenceToTranslate.replace(/(...)/g, function (
          str,
          p1,
          offset,
          s
        ) {
          return " " + p1 + " ";
        });

        for (var m = 0; m < geneticCodeMatchExp.length; m++) {
          sequenceToTranslate = sequenceToTranslate.replace(
            geneticCodeMatchExp[m],
            geneticCodeMatchResult[m]
          );
        }
        sequenceToTranslate = sequenceToTranslate.replace(/\S{3}/g, "X");
        sequenceToTranslate = sequenceToTranslate.replace(/\s/g, "");
        sequenceToTranslate = sequenceToTranslate.replace(/[a-z]/g, "");
        outputWindow.document.write(addReturns(sequenceToTranslate) + "\n\n");

        k = k + 1;
        break;
      }
    }
    startPos = i + 3;
    i = startPos;
    foundStart = false;
    foundStop = false;
    proteinLength = 0;
  }
  if (k == 0) {
    outputWindow.document.write(
      "No ORFs were found in reading frame " + rf + ".\n\n"
    );
  }
  return true;
}
