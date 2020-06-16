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

function emblFasta(theDocument) {
  var newDna = "";
  var maxInput = 200000000;

  if (testScript() == false) {
    return false;
  }

  if (
    checkFormElement(theDocument.forms[0].elements[0]) == false ||
    verifyEmbl(theDocument.forms[0].elements[0].value) == false ||
    checkTextLength(theDocument.forms[0].elements[0].value, maxInput) == false
  ) {
    return false;
  }

  openWindow("EMBL to FASTA");
  openPre();
  emblToFasta(theDocument.forms[0].elements[0].value);
  closePre();
  closeWindow();
  return true;
}

function emblToFasta(emblFile) {
  var title;
  emblFile = "_" + emblFile + "_";
  var recordArray = emblFile.split(/ID\s\s\s[^\f\n\r]*/);
  for (var i = 1; i < recordArray.length; i++) {
    var mainArray = recordArray[i].split(
      /[\f\n\r]\s*FH   Key             Location\/Qualifiers[\f\n\r]+\s*FH|[\f\n\r]\s*XX[\s]*[\f\n\r]\s*SQ[^\f\n\r]*/
    );
    if (mainArray[0].search(/[\f\n\r]\s*DE[^\f\n\r]+/) != -1) {
      title = mainArray[0]
        .match(/[\f\n\r]\s*DE[^\f\n\r]+/)
        .toString()
        .replace(/[\f\n\r]\s*DE\s*/, "");
    } else {
      title = "Untitled";
    }
    title = filterFastaTitle(title.replace(/[\f\n\r\t]+$/g, "")) + "\n";
    dnaArray = mainArray[2].split(/\/{2}/);
    if (dnaArray.length == 1) {
      alert("The entire EMBL file may not have been processed.");
      outputWindow.focus();
    }
    dnaSequence = removeNonDna(dnaArray[0]);
    outputWindow.document.write(
      "&gt;" + title + addReturns(dnaSequence) + "\n\n"
    );
  }
  return true;
}
