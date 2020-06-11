//Written by Paul Stothard, University of Alberta, Canada

function shuffleProtein(theDocument) {
  var newProtein = "";
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

  openWindow("Shuffle Protein");
  openPre();
  var arrayOfFasta = getArrayOfFasta(theDocument.forms[0].elements[0].value);

  for (var i = 0; i < arrayOfFasta.length; i++) {
    newProtein = getSequenceFromFasta(arrayOfFasta[i]);
    title = getTitleFromFasta(arrayOfFasta[i]);

    newProtein = removeNonProteinAllowDegen(newProtein);

    outputWindow.document.write(
      getFastaTitleFromTitleAndSequence(title, newProtein)
    );

    writeShuffledSequence(newProtein);

    outputWindow.document.write("\n\n");
  }

  closePre();
  closeWindow();
  return true;
}
