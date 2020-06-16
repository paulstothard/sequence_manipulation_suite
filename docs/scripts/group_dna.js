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

function groupDna(theDocument) {
  var newDna = "";
  var title = "";
  var maxInput = 100000000;
  var adjustedStart;

  if (testScript() == false) {
    return false;
  }

  if (
    checkFormElement(theDocument.forms[0].elements[0]) == false ||
    checkSequenceLength(theDocument.forms[0].elements[0].value, maxInput) ==
      false
  ) {
    return false;
  }

  adjustedStart = theDocument.forms[0].elements[8].value.replace(
    /[^\d\-]/g,
    ""
  );
  if (
    checkFormElement(theDocument.forms[0].elements[8]) == false ||
    verifyMaxDigits(adjustedStart, 9999999999) == false
  ) {
    return false;
  }
  adjustedStart = parseInt(adjustedStart) - 1;

  openWindow("Group DNA");
  openPre();
  var arrayOfFasta = getArrayOfFasta(theDocument.forms[0].elements[0].value);

  for (var i = 0; i < arrayOfFasta.length; i++) {
    newDna = getSequenceFromFasta(arrayOfFasta[i]);
    title = getTitleFromFasta(arrayOfFasta[i]);

    newDna = removeNonDna(newDna);

    outputWindow.document.write(getInfoFromTitleAndSequence(title, newDna));

    writeGroupNumDnaSetStart(
      newDna,
      "",
      theDocument.forms[0].elements[4].options[
        theDocument.forms[0].elements[4].selectedIndex
      ].value,
      theDocument.forms[0].elements[5].options[
        theDocument.forms[0].elements[5].selectedIndex
      ].value,
      0,
      newDna.length,
      theDocument.forms[0].elements[6].options[
        theDocument.forms[0].elements[6].selectedIndex
      ].value,
      theDocument.forms[0].elements[7].options[
        theDocument.forms[0].elements[7].selectedIndex
      ].value,
      adjustedStart
    );
    outputWindow.document.write("\n\n");
  }

  closePre();
  closeWindow();
  return true;
}
