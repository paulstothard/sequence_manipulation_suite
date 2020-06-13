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

function genbankFasta(theDocument) {
  var newDna = "";
  var maxInput = 200000000;

  if (testScript() == false) {
    return false;
  }

  if (
    checkFormElement(theDocument.forms[0].elements[0]) == false ||
    verifyGenBank(theDocument.forms[0].elements[0].value) == false ||
    checkTextLength(theDocument.forms[0].elements[0].value, maxInput) == false
  ) {
    return false;
  }

  openWindow("GenBank to FASTA");
  openPre();
  genbankToFasta(theDocument.forms[0].elements[0].value);
  closePre();
  closeWindow();
  return true;
}

function genbankToFasta(genBankFile) {
  genBankFile = "_" + genBankFile + "_";
  var recordArray = genBankFile.split(/LOCUS\s\s\s[^\f\n\r]*/m);
  for (var i = 1; i < recordArray.length; i++) {
    var mainArray = recordArray[i].split(
      /DEFINITION|ACCESSION|ORIGIN[^\f\n\r]*/
    );
    var title =
      filterFastaTitle(mainArray[1].replace(/[\f\n\r\t]+$/g, "")) + "\n";
    var dnaSequenceArray = mainArray[3].split(/\/{2}/);
    if (dnaSequenceArray.length == 1) {
      alert("The entire GenBank file may not have been processed.");
      outputWindow.focus();
    }
    var dnaSequence = removeNonDna(dnaSequenceArray[0]);
    outputWindow.document.write(
      "&gt;" + title + addReturns(dnaSequence) + "\n\n"
    );
  }
  return true;
}
