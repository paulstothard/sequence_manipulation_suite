//Written by Paul Stothard, University of Alberta, Canada

function proteinPattern(theDocument) {
  var newProtein = "";
  var maxInput = 500000000;

  if (testScript() == false) {
    return false;
  }

  var re;
  var matches = new Array();

  if (
    checkFormElement(theDocument.forms[0].elements[0]) == false ||
    checkSequenceLength(theDocument.forms[0].elements[0].value, maxInput) ==
      false ||
    checkFormElement(theDocument.forms[0].elements[1]) == false
  ) {
    return false;
  }

  var re =
    "/" + theDocument.forms[0].elements[1].value.replace(/\//g, "") + "/gi";
  re = removeWhiteSpace(re);
  try {
    re = eval(re);
    var testString = "teststring";
    testString = testString.replace(re, "");
  } catch (e) {
    alert("The regular expression is not formatted correctly.");
    return false;
  }

  openWindow("Protein Pattern Find");
  openPre();
  var arrayOfFasta = getArrayOfFasta(theDocument.forms[0].elements[0].value);

  for (var i = 0; i < arrayOfFasta.length; i++) {
    newProtein = getSequenceFromFasta(arrayOfFasta[i]);
    title = getTitleFromFasta(arrayOfFasta[i]);
    newProtein = removeNonProteinStrict(newProtein);

    outputWindow.document.write(getInfoFromTitleAndSequence(title, newProtein));

    writeProteinPattern(newProtein, re);

    outputWindow.document.write("\n\n");
  }

  closePre();
  closeWindow();
  return true;
}

function writeProteinPattern(proteinSequence, re) {
  var matchArray;
  var matchPosition;
  var matchCount = 0;
  var simplePattern = re.toString();
  simplePattern = simplePattern.replace(/\/gi$|\/ig$|^\//gi, "");

  while ((matchArray = re.exec(proteinSequence))) {
    matchCount++;

    var match_end = re.lastIndex;
    var match_start = match_end - RegExp.lastMatch.length + 1;

    outputWindow.document.write(
      "&gt;match number " +
        matchCount +
        ' to "' +
        simplePattern +
        '" start=' +
        match_start +
        " end=" +
        match_end +
        "\n" +
        addReturns(matchArray[0]) +
        "\n\n"
    );

    re.lastIndex = re.lastIndex - RegExp.lastMatch.length + 1;
  }
  if (!(matchCount > 0)) {
    outputWindow.document.write("no matches found for this sequence.\n\n");
  }
}
