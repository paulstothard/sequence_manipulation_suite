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

function dnaPattern(theDocument) {
  var newDna = "";
  var maxInput = 500000000;
  var matches = new Array();

  if (testScript() == false) {
    return false;
  }

  if (
    checkFormElement(theDocument.forms[0].elements[0]) == false ||
    checkSequenceLength(theDocument.forms[0].elements[0].value, maxInput) ==
      false ||
    checkFormElement(theDocument.forms[0].elements[1]) == false
  ) {
    return false;
  }

  var re =
    "/" + theDocument.forms[0].elements[1].value.replace(/\//g, "") + "/gi";
  re = removeWhiteSpace(re);
  try {
    re = eval(re);
    var testString = "teststring";
    testString = testString.replace(re, "");
  } catch (e) {
    alert("The regular expression is not formatted correctly.");
    return false;
  }

  openWindow("DNA Pattern Find");
  openPre();
  var arrayOfFasta = getArrayOfFasta(theDocument.forms[0].elements[0].value);

  for (var i = 0; i < arrayOfFasta.length; i++) {
    newDna = getSequenceFromFasta(arrayOfFasta[i]);
    title = getTitleFromFasta(arrayOfFasta[i]);
    newDna = removeNonDna(newDna);

    outputWindow.document.write(getInfoFromTitleAndSequence(title, newDna));

    writeDnaPattern(newDna, re);

    outputWindow.document.write("\n\n");
  }

  closePre();
  closeWindow();
  return true;
}

function writeDnaPattern(dnaSequence, re) {
  var matchArray;
  var matchCount = 0;
  var length = dnaSequence.length;
  var simplePattern = re.toString();
  simplePattern = simplePattern.replace(/\/gi$|\/ig$|^\//gi, "");

  while ((matchArray = re.exec(dnaSequence))) {
    matchCount++;
    var match_end = re.lastIndex;
    var match_start = match_end - RegExp.lastMatch.length + 1;

    outputWindow.document.write(
      "&gt;match number " +
        matchCount +
        ' to "' +
        simplePattern +
        '" start=' +
        match_start +
        " end=" +
        match_end +
        " on the direct strand\n" +
        addReturns(matchArray[0]) +
        "\n\n"
    );

    re.lastIndex = re.lastIndex - RegExp.lastMatch.length + 1;
  }

  //reset search to start of sequence
  re.lastIndex = 0;
  //now search the reverse-complement
  dnaSequence = reverse(complement(dnaSequence));
  while ((matchArray = re.exec(dnaSequence))) {
    matchCount++;
    var match_start = length - re.lastIndex + 1;
    var match_end = match_start + RegExp.lastMatch.length - 1;

    outputWindow.document.write(
      "&gt;match number " +
        matchCount +
        ' to "' +
        simplePattern +
        '" start=' +
        match_start +
        " end=" +
        match_end +
        " on the reverse strand\n" +
        addReturns(matchArray[0]) +
        "\n\n"
    );
    re.lastIndex = re.lastIndex - RegExp.lastMatch.length + 1;
  }

  if (!(matchCount > 0)) {
    outputWindow.document.write("no matches found for this sequence.\n\n");
  }
}
