//Written by Paul Stothard, University of Alberta, Canada

function shuffleDna(theDocument) {
  var newDna = "";
  var title = "";
  var maxInput = 300000000;

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

  openWindow("Shuffle DNA");
  openPre();
  var arrayOfFasta = getArrayOfFasta(theDocument.forms[0].elements[0].value);

  for (var i = 0; i < arrayOfFasta.length; i++) {
    newDna = getSequenceFromFasta(arrayOfFasta[i]);
    title = getTitleFromFasta(arrayOfFasta[i]);

    newDna = removeNonDna(newDna);

    outputWindow.document.write(
      getFastaTitleFromTitleAndSequence(title, newDna)
    );

    writeShuffledSequence(newDna);

    outputWindow.document.write("\n\n");
  }

  closePre();
  closeWindow();
  return true;
}
