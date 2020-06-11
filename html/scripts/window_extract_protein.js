//Written by Paul Stothard, University of Alberta, Canada

function windowExtract(theDocument) {
  var newDna = "";
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
    checkFormElement(theDocument.forms[0].elements[1]) == false ||
    checkFormElement(theDocument.forms[0].elements[3]) == false ||
    verifyMaxDigits(
      theDocument.forms[0].elements[1].value.replace(/[^\d]/g, ""),
      maxInput
    ) == false ||
    verifyMaxDigits(
      theDocument.forms[0].elements[3].value.replace(/[^\d]/g, ""),
      maxInput
    ) == false
  ) {
    return false;
  }

  //build single range
  var windowSize = parseInt(
    theDocument.forms[0].elements[1].value.replace(/[^\d]/g, "")
  );
  var position = parseInt(
    theDocument.forms[0].elements[3].value.replace(/[^\d]/g, "")
  );
  var orientation = theDocument.forms[0].elements[2].value;

  var start;
  var end;
  if (orientation == "ending") {
    end = position;
    start = end - windowSize + 1;
  } else if (orientation == "starting") {
    start = position;
    end = start + windowSize - 1;
  } else if (orientation == "centered") {
    start = position - Math.round(windowSize / 2) + 1;
    end = start + windowSize - 1;
  }

  ranges.push(new Range(start, end));

  openWindow("Window Extractor Protein");
  openPre();
  var arrayOfFasta = getArrayOfFasta(theDocument.forms[0].elements[0].value);
  for (var i = 0; i < arrayOfFasta.length; i++) {
    newProtein = getSequenceFromFasta(arrayOfFasta[i]);
    title = getTitleFromFasta(arrayOfFasta[i]);
    newProtein = removeNonProteinAllowDegen(newProtein);
    outputWindow.document.write(
      getFastaTitleFromTitleAndSequence(title, newProtein)
    );
    writeSequenceRanges(
      newProtein,
      ranges,
      theDocument.forms[0].elements[7].options[
        theDocument.forms[0].elements[7].selectedIndex
      ].value
    );
  }
  closePre();
  closeWindow();
  return true;
}
