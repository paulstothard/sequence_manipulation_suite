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
