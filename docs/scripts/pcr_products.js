//    Sequence Manipulation Suite. A collection of simple JavaScript programs
//    for generating, formatting, and analyzing short DNA and protein
//    sequences.
//    Copyright (C) 2020 Paul Stothard stothard@ualberta.ca
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <https://www.gnu.org/licenses/>.
//

//Written by Paul Stothard, University of Alberta, Canada

function pcrProducts(theDocument) {
  var newDna = "";
  var maxInput = 200000000;

  if (testScript() == false) {
    return false;
  }

  var re;
  var primers = new Array();
  var forwardMatches = new Array();
  var reverseMatches = new Array();
  var pcrProducts = new Array();

  if (
    checkFormElement(theDocument.forms[0].elements[0]) == false ||
    checkSequenceLength(theDocument.forms[0].elements[0].value, maxInput) ==
      false ||
    checkFormElement(theDocument.forms[0].elements[1]) == false ||
    checkFormElement(theDocument.forms[0].elements[2]) == false ||
    checkFormElement(theDocument.forms[0].elements[3]) == false ||
    checkFormElement(theDocument.forms[0].elements[4]) == false
  ) {
    return false;
  }

  if (
    theDocument.forms[0].elements[2].value.replace(/[^A-Za-z]/g, "").length < 10
  ) {
    alert("Please enter primer sequences that are at least 10 bases long.");
    return false;
  }

  if (
    theDocument.forms[0].elements[4].value.replace(/[^A-Za-z]/g, "").length < 10
  ) {
    alert("Please enter primer sequences that are at least 10 bases long.");
    return false;
  }

  var primerOne = convertDegenerates(
    theDocument.forms[0].elements[2].value.replace(/[^A-Za-z]/g, "")
  );
  var primerTwo = convertDegenerates(
    theDocument.forms[0].elements[4].value.replace(/[^A-Za-z]/g, "")
  );

  try {
    re = eval("/" + primerOne + "/gi");
    var testString = "teststring";
    testString = testString.replace(re, "");
  } catch (e) {
    alert("The first primer sequence is not formatted correctly.");
    return false;
  }

  try {
    re = eval("/" + primerTwo + "/gi");
    var testString = "teststring";
    testString = testString.replace(re, "");
  } catch (e) {
    alert("The second primer sequence is not formatted correctly.");
    return false;
  }

  primers.push(
    new Primer(
      eval("/" + primerOne + "/gi"),
      theDocument.forms[0].elements[1].value
    )
  );
  primers.push(
    new Primer(
      eval("/" + primerTwo + "/gi"),
      theDocument.forms[0].elements[3].value
    )
  );

  var arrayOfFasta = getArrayOfFasta(theDocument.forms[0].elements[0].value);
  for (var i = 0; i < arrayOfFasta.length; i++) {
    newDna = getSequenceFromFasta(arrayOfFasta[i]);
    title = getTitleFromFasta(arrayOfFasta[i]);
    newDna = removeNonDna(newDna);

    //create matchedPrimer objects from sequences and add them to array
    forwardMatches = findMatches(
      primers,
      newDna,
      theDocument.forms[0].elements[8].options[
        theDocument.forms[0].elements[8].selectedIndex
      ].value
    );
    reverseMatches = findMatches(
      primers,
      reverse(complement(newDna)),
      theDocument.forms[0].elements[8].options[
        theDocument.forms[0].elements[8].selectedIndex
      ].value
    );
    makePcrProducts(
      newDna,
      title,
      forwardMatches,
      reverseMatches,
      theDocument.forms[0].elements[8].options[
        theDocument.forms[0].elements[8].selectedIndex
      ].value,
      pcrProducts
    );
  }

  pcrProducts.sort(pcrProductSorter);

  openWindow("PCR Products");
  openPre();
  for (var i = 0; i < pcrProducts.length; i++) {
    pcrProducts[i].writeProduct(
      theDocument.forms[0].elements[8].options[
        theDocument.forms[0].elements[8].selectedIndex
      ].value
    );
  }
  if (pcrProducts.length == 0) {
    outputWindow.document.write("No PCR products were obtained.\n\n");
  }
  closePre();
  closeWindow();
  return true;
}

function findMatches(primers, sequence, topology) {
  var matchArray;
  var matchPosition;
  var arrayOfMatches = new Array();
  var re;
  var originalLength = sequence.length;

  if (topology == "circular") {
    var lookAhead = 50;
    var shiftValue = sequence.substring(0, lookAhead).length;
    var upperLimit = sequence.length + shiftValue;
    sequence =
      sequence.substring(sequence.length - lookAhead, sequence.length) +
      sequence +
      sequence.substring(0, lookAhead);
    var lowerLimit = 0 + shiftValue;

    for (var i = 0; i < primers.length; i++) {
      re = primers[i].re;
      while ((matchArray = re.exec(sequence))) {
        matchPosition = re.lastIndex;
        if (matchPosition >= lowerLimit && matchPosition < upperLimit) {
          matchPosition = matchPosition - shiftValue;
          if (matchPosition == 0) {
            matchPosition = originalLength;
          }
          arrayOfMatches.push(
            new Match(primers[i].name, matchArray[0], matchPosition)
          );
        }
        re.lastIndex = re.lastIndex - RegExp.lastMatch.length + 1;
      }
    }
  } else {
    for (var i = 0; i < primers.length; i++) {
      re = primers[i].re;
      while ((matchArray = re.exec(sequence))) {
        matchPosition = re.lastIndex;
        arrayOfMatches.push(
          new Match(primers[i].name, matchArray[0], matchPosition)
        );
        re.lastIndex = re.lastIndex - RegExp.lastMatch.length + 1;
      }
    }
  }
  return arrayOfMatches;
}

function makePcrProducts(
  newDna,
  title,
  forwardMatches,
  reverseMatches,
  topology,
  pcrProducts
) {
  for (var i = 0; i < forwardMatches.length; i++) {
    for (var j = 0; j < reverseMatches.length; j++) {
      if (
        forwardMatches[i].positionAfter -
          forwardMatches[i].matchingText.length <=
        newDna.length - reverseMatches[j].positionAfter
      ) {
        pcrProducts.push(
          new PcrProduct(
            title,
            forwardMatches[i].positionAfter -
              forwardMatches[i].matchingText.length +
              1,
            newDna.length -
              reverseMatches[j].positionAfter +
              reverseMatches[j].matchingText.length,
            forwardMatches[i].primerName,
            reverseMatches[j].primerName,
            forwardMatches[i].name,
            reverseMatches[j].name,
            newDna.substring(
              forwardMatches[i].positionAfter -
                forwardMatches[i].matchingText.length,
              newDna.length -
                reverseMatches[j].positionAfter +
                reverseMatches[j].matchingText.length
            )
          )
        );
      } else if (
        topology == "circular" &&
        forwardMatches[i].positionAfter -
          forwardMatches[i].matchingText.length >
          newDna.length -
            reverseMatches[j].positionAfter +
            reverseMatches[j].matchingText.length -
            1
      ) {
        pcrProducts.push(
          new PcrProduct(
            title,
            forwardMatches[i].positionAfter -
              forwardMatches[i].matchingText.length +
              1,
            newDna.length -
              reverseMatches[j].positionAfter +
              reverseMatches[j].matchingText.length,
            forwardMatches[i].primerName,
            reverseMatches[j].primerName,
            forwardMatches[i].name,
            reverseMatches[j].name,
            newDna.substring(
              forwardMatches[i].positionAfter -
                forwardMatches[i].matchingText.length,
              newDna.length
            ) +
              newDna.substring(
                0,
                newDna.length -
                  reverseMatches[j].positionAfter +
                  reverseMatches[j].matchingText.length
              )
          )
        );
      }
    }
  }
}

//PcrProduct class writeProduct() method
function writeProduct(topology) {
  outputWindow.document.write(
    "&gt;" +
      this.sequence.length +
      " bp product from " +
      topology +
      " template " +
      this.template +
      ", base " +
      this.start +
      " to base " +
      this.stop +
      " (" +
      this.forwardName +
      " - " +
      this.reverseName +
      ").\n"
  );
  outputWindow.document.write(addReturns(this.sequence) + "\n\n");
}

//PcrProduct class
function PcrProduct(
  template,
  start,
  stop,
  forwardName,
  reverseName,
  forwardPrimer,
  reversePrimer,
  sequence
) {
  this.template = template;
  this.start = start;
  this.stop = stop;
  this.forwardName = forwardName;
  this.reverseName = reverseName;
  this.forwardPrimer = forwardPrimer;
  this.reversePrimer = reversePrimer;
  this.sequence = sequence;
}

//create and throw away a prototype object
new PcrProduct("", 0, 0, "", "", "", "", "");

// define object methods
PcrProduct.prototype.writeProduct = writeProduct;

//Match class
function Match(primerName, matchingText, positionAfter) {
  this.primerName = primerName;
  this.matchingText = matchingText;
  this.positionAfter = positionAfter;
}

//Primer class
function Primer(re, name) {
  this.re = re;
  this.name = name;
}

function pcrProductSorter(a, b) {
  if (a.sequence.length < b.sequence.length) {
    return 1;
  }
  if (a.sequence.length > b.sequence.length) {
    return -1;
  } else {
    return 0;
  }
}
