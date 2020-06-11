//Written by Paul Stothard, University of Alberta, Canada

function filterDna(theDocument) {
  var newDna = "";
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
  newDna = theDocument.forms[0].elements[0].value.replace(
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
    newDna = newDna.toUpperCase();
  } else {
    if (
      theDocument.forms[0].elements[6].options[
        theDocument.forms[0].elements[6].selectedIndex
      ].value == "lowercase"
    ) {
      newDna = newDna.toLowerCase();
    }
  }

  openWindow("Filter DNA");
  openPre();
  outputWindow.document.write(
    "&gt;filtered DNA sequence consisting of " + newDna.length + " bases.\n"
  );
  outputWindow.document.write(addReturns(newDna));
  outputWindow.document.write("\n");
  closePre();
  closeWindow();
  return true;
}
