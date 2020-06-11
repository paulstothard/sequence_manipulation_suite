//Written by Paul Stothard, University of Alberta, Canada

function sampleProtein(theDocument) {
  var maxInput = 100000000;

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

  var newProtein = getSequenceFromFasta(theDocument.forms[0].elements[1].value);
  var title = getTitleFromFasta(theDocument.forms[0].elements[1].value);
  newProtein = removeNonProteinAllowDegen(newProtein);

  var components = new Array(newProtein.length);
  if (newProtein.search(/./) != -1) {
    components = newProtein.match(/./g);
  }

  openWindow("Sample Protein");
  openPre();
  for (var i = 1; i <= seqNum; i++) {
    outputWindow.document.write(
      "&gt;" +
        "sampled sequence " +
        i +
        " consisting of " +
        enteredNumber +
        " residues.\n"
    );
    writeRandomSequence(components, enteredNumber);
    outputWindow.document.write("\n");
  }
  closePre();
  closeWindow();
  return true;
}
