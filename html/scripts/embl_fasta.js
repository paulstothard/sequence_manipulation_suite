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
