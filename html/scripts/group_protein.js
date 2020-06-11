//Written by Paul Stothard, University of Alberta, Canada

function groupProtein(theDocument) {
  var newProtein = "";
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

  openWindow("Group Protein");
  openPre();
  var arrayOfFasta = getArrayOfFasta(theDocument.forms[0].elements[0].value);

  for (var i = 0; i < arrayOfFasta.length; i++) {
    newProtein = getSequenceFromFasta(arrayOfFasta[i]);
    title = getTitleFromFasta(arrayOfFasta[i]);

    newProtein = removeNonProtein(newProtein);

    outputWindow.document.write(getInfoFromTitleAndSequence(title, newProtein));

    writeGroupNumProtein(
      newProtein,
      "",
      theDocument.forms[0].elements[4].options[
        theDocument.forms[0].elements[4].selectedIndex
      ].value,
      theDocument.forms[0].elements[5].options[
        theDocument.forms[0].elements[5].selectedIndex
      ].value,
      0,
      newProtein.length,
      theDocument.forms[0].elements[6].options[
        theDocument.forms[0].elements[6].selectedIndex
      ].value
    );
    outputWindow.document.write("\n\n");
  }

  closePre();
  closeWindow();
  return true;
}
