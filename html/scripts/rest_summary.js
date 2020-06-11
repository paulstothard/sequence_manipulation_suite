//Written by Paul Stothard, University of Alberta, Canada

function restSummary(theDocument) {
  var newDna = "";
  var title = "";
  var maxInput = 100000000;

  var restrictionSites = getRestrictionSiteString("standard");

  if (
    checkFormElement(theDocument.forms[0].elements[0]) == false ||
    checkSequenceLength(theDocument.forms[0].elements[0].value, maxInput) ==
      false
  ) {
    return false;
  }

  itemsToCheck = restrictionSites.split(/,/);
  if (checkRestPatterns(itemsToCheck) == false) {
    return false;
  }

  openWindow("Restriction Summary");
  //openPre();
  outputWindow.document.write(
    '<span class="one">' + "cuts once" + "</span><br />\n"
  );
  outputWindow.document.write(
    '<span class="two">' + "cuts twice" + "</span><br />\n"
  );
  outputWindow.document.write("\n");
  var arrayOfFasta = getArrayOfFasta(theDocument.forms[0].elements[0].value);

  for (var i = 0; i < arrayOfFasta.length; i++) {
    newDna = getSequenceFromFasta(arrayOfFasta[i]);
    title = getTitleFromFasta(arrayOfFasta[i]);

    newDna = removeNonDna(newDna);

    outputWindow.document.write(
      getInfoFromTitleAndSequenceAndTopology(
        title,
        newDna,
        theDocument.forms[0].elements[4].options[
          theDocument.forms[0].elements[4].selectedIndex
        ].value
      )
    );

    writeRestrictionSites(
      newDna,
      itemsToCheck,
      theDocument.forms[0].elements[4].options[
        theDocument.forms[0].elements[4].selectedIndex
      ].value
    );

    outputWindow.document.write("<br />\n<br />\n");
  }

  //closePre();
  closeWindow();
  return true;
}
