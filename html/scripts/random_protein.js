//Written by Paul Stothard, University of Alberta, Canada

function randomProtein(theDocument) {
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
  openWindow("Random Protein Sequence");
  openPre();
  for (var i = 1; i <= seqNum; i++) {
    outputWindow.document.write(
      "&gt;" +
        "random sequence " +
        i +
        " consisting of " +
        enteredNumber +
        " residues.\n"
    );
    writeRandomSequence(
      [
        "A",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y",
      ],
      enteredNumber
    );
    outputWindow.document.write("\n");
  }
  closePre();
  closeWindow();
  return true;
}
