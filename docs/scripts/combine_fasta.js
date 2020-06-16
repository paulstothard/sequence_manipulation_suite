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

function combineFasta(theDocument) {
  var maxInput = 500000000;
  var sequenceCount = 0;
  var sequences = new Array();

  if (testScript() == false) {
    return false;
  }

  if (
    checkFormElement(theDocument.forms[0].elements[0]) == false ||
    checkTextLength(theDocument.forms[0].elements[0].value, maxInput) == false
  ) {
    return false;
  }

  var arrayOfFasta = getArrayOfFasta(theDocument.forms[0].elements[0].value);

  for (var i = 0; i < arrayOfFasta.length; i++) {
    sequences.push(removeNonLetters(getSequenceFromFasta(arrayOfFasta[i])));
  }

  var sequence = sequences.join("");
  openWindow("Combine FASTA");
  openPre();
  if (sequences.length == 1) {
    outputWindow.document.write(
      "&gt;results for " +
        sequence.length +
        " residue sequence made from " +
        sequences.length +
        ' records, starting "' +
        sequence.substring(0, 10) +
        '"\n'
    );
  } else if (sequences.length > 1) {
    outputWindow.document.write(
      "&gt;results for " +
        sequence.length +
        " residue sequence made from " +
        sequences.length +
        ' records, starting "' +
        sequence.substring(0, 10) +
        '"\n'
    );
  } else {
    outputWindow.document.write(
      '<div class="info">No sequence records were read</div>\n'
    );
  }
  outputWindow.document.write(addReturns(sequence) + "\n");
  closePre();
  closeWindow();
  return true;
}
