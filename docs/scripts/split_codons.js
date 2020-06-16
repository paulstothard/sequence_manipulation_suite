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

function splitCodons(theDocument) {
  var maxInput = 500000000;
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

  openWindow("Split Codons");
  openPre();
  var arrayOfFasta = getArrayOfFasta(theDocument.forms[0].elements[0].value);

  for (var i = 0; i < arrayOfFasta.length; i++) {
    var sequence = getSequenceFromFasta(arrayOfFasta[i]);
    sequence = removeFormatting(sequence);
    var title = getTitleFromFasta(arrayOfFasta[i]);

    if (sequence.length % 3 != 0) {
      alert(
        "Sequence '" + title + "' ends in a partial codon that will be removed."
      );
    }

    var length = sequence.length;
    var seqCount = 1;

    var position1 = getBasesBasedOnCodonPosition(sequence, 1);
    outputWindow.document.write(
      ">" +
        title +
        ";codon_positon_1_bases;length=" +
        position1.length +
        ";source_length=" +
        length +
        "\n" +
        addReturns(position1) +
        "\n\n"
    );

    var position2 = getBasesBasedOnCodonPosition(sequence, 2);
    outputWindow.document.write(
      ">" +
        title +
        ";codon_positon_2_bases;length=" +
        position2.length +
        ";source_length=" +
        length +
        "\n" +
        addReturns(position2) +
        "\n\n"
    );

    var position3 = getBasesBasedOnCodonPosition(sequence, 3);
    outputWindow.document.write(
      ">" +
        title +
        ";codon_positon_3_bases;length=" +
        position3.length +
        ";source_length=" +
        length +
        "\n" +
        addReturns(position3) +
        "\n\n"
    );

    seqCount++;
  }

  closePre();
  closeWindow();
  return true;
}

function getBasesBasedOnCodonPosition(sequence, position) {
  var re;
  if (position == 1) {
    re = "((.)..)";
  } else if (position == 2) {
    re = "(.(.).)";
  } else if (position == 3) {
    re = "(..(.))";
  }

  //remove partial codon from the end
  var partial_codon_length = sequence.length % 3;
  sequence = sequence.replace(
    new RegExp(".{" + partial_codon_length + "}$"),
    ""
  );

  return sequence.replace(new RegExp(re, "g"), function (
    str,
    p1,
    p2,
    offset,
    s
  ) {
    return p2;
  });
}
