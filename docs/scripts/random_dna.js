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

function randomDna(theDocument) {
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
    verifyMaxDigits(enteredNumber, maxInput) == false
  ) {
    return false;
  }

  var seqNum = parseInt(
    theDocument.forms[0].elements[4].options[
      theDocument.forms[0].elements[4].selectedIndex
    ].value
  );
  openWindow("Random DNA Sequence");
  openPre();
  for (var i = 1; i <= seqNum; i++) {
    outputWindow.document.write(
      "&gt;" +
        "random sequence " +
        i +
        " consisting of " +
        enteredNumber +
        " bases.\n"
    );
    writeRandomSequence(["g", "a", "c", "t"], enteredNumber);
    outputWindow.document.write("\n");
  }
  closePre();
  closeWindow();
  return true;
}
