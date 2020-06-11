//Written by Paul Stothard, University of Alberta, Canada

function mutateDna(theDocument) {
  var newDna = "";
  var title = "";
  var maxInput = 1000000000;
  var maxDigitsInput = 100000000;

  if (testScript() == false) {
    return false;
  }

  var enteredNumber = theDocument.forms[0].elements[4].value.replace(
    /[^\d]/g,
    ""
  );

  if (
    checkFormElement(theDocument.forms[0].elements[0]) == false ||
    checkSequenceLength(theDocument.forms[0].elements[0].value, maxInput) ==
      false ||
    verifyMaxDigits(enteredNumber, maxDigitsInput) == false
  ) {
    return false;
  }

  openWindow("Mutate DNA");
  openPre();
  var arrayOfFasta = getArrayOfFasta(theDocument.forms[0].elements[0].value);

  for (var i = 0; i < arrayOfFasta.length; i++) {
    newDna = getSequenceFromFasta(arrayOfFasta[i]);
    title = getTitleFromFasta(arrayOfFasta[i]);

    newDna = removeNonDna(newDna);

    outputWindow.document.write(
      getFastaTitleFromTitleAndSequence(title, newDna)
    );
    writeMutatedSequence(
      newDna,
      ["g", "a", "c", "t"],
      enteredNumber,
      theDocument.forms[0].elements[5].options[
        theDocument.forms[0].elements[5].selectedIndex
      ].value,
      newDna.length -
        theDocument.forms[0].elements[6].options[
          theDocument.forms[0].elements[6].selectedIndex
        ].value -
        1
    );

    outputWindow.document.write("\n\n");
  }

  closePre();
  closeWindow();
  return true;
}
