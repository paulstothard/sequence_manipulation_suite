//Written by Paul Stothard, University of Alberta, Canada

function randomDnaRegions(theDocument) {
  var newDna = "";
  var title = "";
  var maxInput = 500000000;
  var matchFound = false;
  var ranges = new Array();

  if (testScript() == false) {
    return false;
  }

  if (
    checkFormElement(theDocument.forms[0].elements[0]) == false ||
    checkSequenceLength(theDocument.forms[0].elements[0].value, maxInput) ==
      false ||
    checkFormElement(theDocument.forms[0].elements[1]) == false
  ) {
    return false;
  }

  var arrayOfRanges = theDocument.forms[0].elements[1].value.split(/,/);
  var arrayOfStartAndEnd;
  for (var i = 0; i < arrayOfRanges.length; i++) {
    arrayOfStartAndEnd = arrayOfRanges[i].split(/\.\.|\-/);
    if (
      arrayOfStartAndEnd.length == 1 &&
      arrayOfStartAndEnd[0].search(/\d/) != -1
    ) {
      matchFound = true;
      ranges.push(new Range(arrayOfStartAndEnd[0], arrayOfStartAndEnd[0]));
    } else if (
      arrayOfStartAndEnd.length == 2 &&
      arrayOfStartAndEnd[0].search(/\d/) != -1 &&
      arrayOfStartAndEnd[1].search(/\d/) != -1
    ) {
      matchFound = true;
      ranges.push(new Range(arrayOfStartAndEnd[0], arrayOfStartAndEnd[1]));
    }
  }
  if (matchFound == false) {
    alert("No ranges were entered.");
    return false;
  }

  openWindow("Random DNA Regions");
  openPre();
  var arrayOfFasta = getArrayOfFasta(theDocument.forms[0].elements[0].value);

  for (var i = 0; i < arrayOfFasta.length; i++) {
    newDna = getSequenceFromFasta(arrayOfFasta[i]);
    title = getTitleFromFasta(arrayOfFasta[i]);

    newDna = removeNonDna(newDna);

    outputWindow.document.write(
      getFastaTitleFromTitleAndSequence(title, newDna)
    );

    writeSequenceRanges(
      newDna,
      ranges,
      "direct",
      theDocument.forms[0].elements[5].options[
        theDocument.forms[0].elements[5].selectedIndex
      ].value
    );

    outputWindow.document.write("\n\n");
  }

  closePre();
  closeWindow();
  return true;
}
