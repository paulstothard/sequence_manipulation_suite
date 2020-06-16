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

function cpgIslands(theDocument) {
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

  openWindow("CpG Islands");
  openPre();
  var arrayOfFasta = getArrayOfFasta(theDocument.forms[0].elements[0].value);

  for (var i = 0; i < arrayOfFasta.length; i++) {
    newDna = getSequenceFromFasta(arrayOfFasta[i]);
    title = getTitleFromFasta(arrayOfFasta[i]);
    newDna = removeNonDna(newDna);
    outputWindow.document.write(getInfoFromTitleAndSequence(title, newDna));
    cpgIslandRegions(newDna, 200, 0.6);
  }

  closePre();
  closeWindow();
  return true;
}

function cpgIslandRegions(dnaSequence, windowSize, cutOff) {
  var islandFound = false;
  var numG = 0;
  var numC = 0;
  var numCG = 0;
  var valueY = 0;
  var gcContent = 0;
  dnaSequence = dnaSequence.toLowerCase();
  windowSize = parseInt(windowSize);
  cutOff = parseFloat(cutOff);

  if (windowSize > dnaSequence.length) {
    outputWindow.document.write(
      "The input sequence must be longer than " + windowSize + " bases.<br />\n"
    );
    return true;
  }

  //determine base counts for first window
  for (var i = 0; i < windowSize; i++) {
    if (dnaSequence.charAt(i) == "g") {
      numG = numG + 1;
    }
    if (dnaSequence.charAt(i) == "c") {
      numC = numC + 1;
      if (dnaSequence.charAt(i + 1) == "g") {
        numCG = numCG + 1;
        numG = numG + 1;
        i = i + 1;
      }
    }
  }
  if (numC != 0 && numG != 0) {
    valueY = numCG / ((numC * numG) / windowSize);
  } else {
    valueY = 0;
  }
  gcContent = (numG + numC) / windowSize;
  if (valueY >= cutOff && gcContent > 0.5) {
    //valueY = Math.round(valueY*100)/100;
    //gcContent = Math.round(gcContent*1000)/10;
    gcContent = gcContent * 100;
    valueY = valueY.toFixed(2);
    gcContent = gcContent.toFixed(2);
    outputWindow.document.write(
      "CpG island detected in region 1 to " +
        windowSize +
        " (Obs/Exp = " +
        valueY +
        " and %GC = " +
        gcContent +
        ")<br />\n"
    );
    islandFound = true;
  }

  //now update values as window slides along at 1 base intervals
  start = windowSize;
  for (var j = start; j < dnaSequence.length; j++) {
    baseToAdd = dnaSequence.charAt(j);
    baseToLose = dnaSequence.charAt(j - windowSize);
    recentBaseAdded = dnaSequence.charAt(j - 1);
    nextToLose = dnaSequence.charAt(j - windowSize + 1);
    if (baseToAdd == "c") {
      numC = numC + 1;
    }
    if (baseToAdd == "g") {
      numG = numG + 1;
      if (recentBaseAdded == "c") {
        numCG = numCG + 1;
      }
    }
    if (baseToLose == "c") {
      numC = numC - 1;
      if (nextToLose == "g") {
        numCG = numCG - 1;
      }
    }
    if (baseToLose == "g") {
      numG = numG - 1;
    }
    if (numC != 0 && numG != 0) {
      valueY = numCG / ((numC * numG) / windowSize);
    } else {
      valueY = 0;
    }
    gcContent = (numG + numC) / windowSize;
    if (valueY > cutOff && gcContent > 0.5) {
      startRange = (j - windowSize + 2).toString();
      endRange = (j + 1).toString();
      //valueY = Math.round(valueY*100)/100;
      //gcContent = Math.round(gcContent*1000)/10;
      gcContent = gcContent * 100;
      valueY = valueY.toFixed(2);
      gcContent = gcContent.toFixed(2);
      outputWindow.document.write(
        "CpG island detected in region " +
          startRange +
          " to " +
          endRange +
          " (Obs/Exp = " +
          valueY +
          " and %GC = " +
          gcContent +
          ") <br />\n"
      );
      islandFound = true;
    }
  }
  if (!islandFound) {
    outputWindow.document.write(
      "No CpG island regions were identified.<br />\n"
    );
  }
  return true;
}
