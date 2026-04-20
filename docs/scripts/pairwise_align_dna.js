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

function pairwiseAlignDna(theDocument) {
  var newDnaOne = "";
  var titleOne = "";
  var newDnaTwo = "";
  var titleTwo = "";
  var maxInput = 20000;
  var gapOpen;
  var gapExtend;
  var alignment;

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

  gapOpen = parseInt(
    theDocument.forms[0].elements[5].options[
      theDocument.forms[0].elements[5].selectedIndex
    ].value,
    10
  );
  gapExtend = parseInt(
    theDocument.forms[0].elements[6].options[
      theDocument.forms[0].elements[6].selectedIndex
    ].value,
    10
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

  alignment = PairwiseAffineAlignment.align(
    newDnaOne,
    newDnaTwo,
    "ednafull",
    gapOpen,
    gapExtend
  );

  outputWindow.document.write(">" + titleOne + "\n");
  outputWindow.document.write(addReturns(alignment.alignedM));
  outputWindow.document.write("\n\n");
  outputWindow.document.write(">" + titleTwo + "\n");
  outputWindow.document.write(addReturns(alignment.alignedN));
  outputWindow.document.write("\n\n");
  outputWindow.document.write("Alignment score: " + alignment.score + "\n\n");

  closePre();
  closeWindow();
  return true;
}
