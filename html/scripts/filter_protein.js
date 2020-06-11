//Written by Paul Stothard, University of Alberta, Canada

function filterProtein(theDocument) {
  var newProtein = "";
  var maxInput = 500000000;

  if (testScript() == false) {
    return false;
  }

  if (
    checkFormElement(theDocument.forms[0].elements[0]) == false ||
    checkTextLength(theDocument.forms[0].elements[0].value, maxInput) == false
  ) {
    return false;
  }

  var re = new RegExp(
    theDocument.forms[0].elements[4].options[
      theDocument.forms[0].elements[4].selectedIndex
    ].value,
    "g"
  );
  newProtein = theDocument.forms[0].elements[0].value.replace(
    re,
    theDocument.forms[0].elements[5].options[
      theDocument.forms[0].elements[5].selectedIndex
    ].value
  );
  if (
    theDocument.forms[0].elements[6].options[
      theDocument.forms[0].elements[6].selectedIndex
    ].value == "uppercase"
  ) {
    newProtein = newProtein.toUpperCase();
  } else {
    if (
      theDocument.forms[0].elements[6].options[
        theDocument.forms[0].elements[6].selectedIndex
      ].value == "lowercase"
    ) {
      newProtein = newProtein.toLowerCase();
    }
  }

  openWindow("Filter Protein");
  openPre();
  outputWindow.document.write(
    "&gt;filtered protein sequence consisting of " +
      newProtein.length +
      " residues.\n"
  );
  outputWindow.document.write(addReturns(newProtein));
  outputWindow.document.write("\n");
  closePre();
  closeWindow();
  return true;
}
