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

function translateDna(theDocument) {
  translate(theDocument);
  return true;
}

function translate(theDocument) {
  var newDna = "";
  var title = "";
  var maxInput = 200000000;

  if (testScript() == false) {
    return false;
  }

  var geneticCode = getGeneticCodeString(
    theDocument.forms[0].elements[6].options[
      theDocument.forms[0].elements[6].selectedIndex
    ].value
  );

  if (
    checkFormElement(theDocument.forms[0].elements[0]) == false ||
    checkSequenceLength(theDocument.forms[0].elements[0].value, maxInput) ==
      false
  ) {
    return false;
  }

  geneticCode = geneticCode.split(/,/);

  if (checkGeneticCode(geneticCode) == false) {
    return false;
  }

  var rfText =
    theDocument.forms[0].elements[4].options[
      theDocument.forms[0].elements[4].selectedIndex
    ].value;
  if (
    theDocument.forms[0].elements[4].options[
      theDocument.forms[0].elements[4].selectedIndex
    ].value.match(/^\d+$/)
  ) {
    rfText++;
  } else {
    rfText =
      '"' +
      theDocument.forms[0].elements[4].options[
        theDocument.forms[0].elements[4].selectedIndex
      ].value +
      '"';
  }

  openWindow("Translate");
  openPre();
  var arrayOfFasta = getArrayOfFasta(theDocument.forms[0].elements[0].value);
  for (var i = 0; i < arrayOfFasta.length; i++) {
    newDna = getSequenceFromFasta(arrayOfFasta[i]);
    title = getTitleFromFasta(arrayOfFasta[i]);
    newDna = removeNonDna(newDna);
    outputWindow.document.write("&gt;rf " + rfText + " " + title + "\n");
    writeTranslation(
      newDna,
      geneticCode,
      theDocument.forms[0].elements[4].options[
        theDocument.forms[0].elements[4].selectedIndex
      ].value,
      theDocument.forms[0].elements[5].options[
        theDocument.forms[0].elements[5].selectedIndex
      ].value
    );
    outputWindow.document.write("\n\n");
  }
  closePre();
  closeWindow();
  return true;
}

function writeTranslation(dnaSequence, geneticCode, startPos, strand) {
  var geneticCodeMatchExp = getGeneticCodeMatchExp(geneticCode);
  var geneticCodeMatchResult = getGeneticCodeMatchResult(geneticCode);

  if (strand == "reverse") {
    dnaSequence = reverse(complement(dnaSequence));
  }
  if (startPos == "uppercase") {
    dnaSequence = dnaSequence.replace(/[a-z]/g, "");
  } else {
    dnaSequence = dnaSequence.substring(parseInt(startPos), dnaSequence.length);
  }

  //don't translate if fewer than three bases
  if (dnaSequence.replace(/[^A-Za-z]/g, "").length < 3) {
    return "";
  }

  dnaSequence = dnaSequence.replace(/(...)/g, function (str, p1, offset, s) {
    return " " + p1 + " ";
  });

  for (var i = 0; i < geneticCodeMatchExp.length; i++) {
    dnaSequence = dnaSequence.replace(
      geneticCodeMatchExp[i],
      geneticCodeMatchResult[i]
    );
  }

  dnaSequence = dnaSequence.replace(/\S{3}/g, "X");
  dnaSequence = dnaSequence.replace(/\s\S{1,2}$/, "");
  dnaSequence = dnaSequence.replace(/\s/g, "");
  outputWindow.document.write(addReturns(dnaSequence));
  return true;
}
