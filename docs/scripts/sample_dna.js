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

function sampleDna(theDocument) {
  var maxInput = 10000000;

  if (testScript() == false) {
    return false;
  }

  var enteredNumber = theDocument.forms[0].elements[0].value.replace(
    /[^\d]/g,
    ""
  );
  if (
    checkFormElement(theDocument.forms[0].elements[0]) == false ||
    verifyMaxDigits(enteredNumber, maxInput) == false ||
    checkFormElement(theDocument.forms[0].elements[1]) == false ||
    checkSequenceLength(theDocument.forms[0].elements[1].value, maxInput) ==
      false
  ) {
    return false;
  }

  var seqNum = parseInt(
    theDocument.forms[0].elements[5].options[
      theDocument.forms[0].elements[5].selectedIndex
    ].value
  );

  var newDna = getSequenceFromFasta(theDocument.forms[0].elements[1].value);
  var title = getTitleFromFasta(theDocument.forms[0].elements[1].value);
  verifyDna(newDna);
  newDna = removeNonDna(newDna);

  var components = new Array(newDna.length);
  if (newDna.search(/./) != -1) {
    components = newDna.match(/./g);
  }

  openWindow("Sample DNA");
  openPre();
  for (var i = 1; i <= seqNum; i++) {
    outputWindow.document.write(
      "&gt;" +
        "sampled sequence " +
        i +
        " consisting of " +
        enteredNumber +
        " bases.\n"
    );
    writeRandomSequence(components, enteredNumber);
    outputWindow.document.write("\n");
  }
  closePre();
  closeWindow();
  return true;
}
