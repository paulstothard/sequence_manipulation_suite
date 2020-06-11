//Written by Paul Stothard, University of Alberta, Canada

function revComp(theDocument) {
  var newDna = "";
  var title = "";
  var maxInput = 100000000;

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

  openWindow("Reverse Complement");
  openPre();
  var arrayOfFasta = getArrayOfFasta(theDocument.forms[0].elements[0].value);

  for (var i = 0; i < arrayOfFasta.length; i++) {
    newDna = getSequenceFromFasta(arrayOfFasta[i]);
    title = getTitleFromFasta(arrayOfFasta[i]);
    newDna = removeNonDna(newDna);

    //outputWindow.document.write(getInfoFromTitleAndSequence(title, newDna));

    if (
      theDocument.forms[0].elements[4].options[
        theDocument.forms[0].elements[4].selectedIndex
      ].value == "reverse-complement"
    ) {
      outputWindow.document.write(">" + title + " reverse complement\n");
      newDna = reverse(complement(newDna));
    } else if (
      theDocument.forms[0].elements[4].options[
        theDocument.forms[0].elements[4].selectedIndex
      ].value == "reverse"
    ) {
      outputWindow.document.write(">" + title + " reverse\n");
      newDna = reverse(newDna);
    } else if (
      theDocument.forms[0].elements[4].options[
        theDocument.forms[0].elements[4].selectedIndex
      ].value == "complement"
    ) {
      outputWindow.document.write(">" + title + " complement\n");
      newDna = complement(newDna);
    }

    outputWindow.document.write(addReturns(newDna) + "\n\n");
  }

  closePre();
  closeWindow();
  return true;
}
