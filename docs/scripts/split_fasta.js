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

function splitFasta(theDocument) {
  var maxInput = 500000000;
  var sequences = new Array();

  if (testScript() == false) {
    return false;
  }

  var newLength = theDocument.forms[0].elements[1].value.replace(/[^\d]/g, "");
  var overlap = theDocument.forms[0].elements[2].value.replace(/[^\d]/g, "");

  if (
    checkFormElement(theDocument.forms[0].elements[0]) == false ||
    checkTextLength(theDocument.forms[0].elements[0].value, maxInput) ==
      false ||
    verifyMaxDigits(newLength, maxInput) == false ||
    verifyMaxDigits(overlap, maxInput) == false
  ) {
    return false;
  }

  newLength = parseInt(newLength);
  overlap = parseInt(overlap);

  openWindow("Split FASTA");
  openPre();
  var arrayOfFasta = getArrayOfFasta(theDocument.forms[0].elements[0].value);

  for (var i = 0; i < arrayOfFasta.length; i++) {
    var sequence = getSequenceFromFasta(arrayOfFasta[i]);
    sequence = removeNonLetters(sequence);
    var title = getTitleFromFasta(arrayOfFasta[i]);
    var length = sequence.length;
    var seqCount = 1;

    for (var j = 0; j < length; j = j + newLength) {
      //if using overlap adjust j
      if (j > overlap) {
        j = j - overlap;
      }
      var subseq = sequence.substring(j, j + newLength);

      var subseq_length = subseq.length;
      var start = j + 1;
      var end = start + subseq_length - 1;
      outputWindow.document.write(
        ">fragment_" +
          seqCount +
          ";" +
          title +
          "_start=" +
          start +
          ";end=" +
          end +
          ";length=" +
          subseq_length +
          ";source_length=" +
          length +
          "\n" +
          addReturns(subseq) +
          "\n\n"
      );
      seqCount++;
    }
  }

  closePre();
  closeWindow();
  return true;
}
