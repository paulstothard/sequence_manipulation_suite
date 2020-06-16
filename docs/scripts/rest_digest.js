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

function restDigest(theDocument) {
  var newDna = "";
  var maxInput = 100000000;

  if (testScript() == false) {
    return false;
  }

  var restrictionFragments = new Array();
  var restrictionFragment;

  if (
    checkFormElement(theDocument.forms[0].elements[0]) == false ||
    checkSequenceLength(theDocument.forms[0].elements[0].value, maxInput) ==
      false
  ) {
    return false;
  }

  var arrayOfFasta = getArrayOfFasta(theDocument.forms[0].elements[0].value);
  for (var i = 0; i < arrayOfFasta.length; i++) {
    newDna = getSequenceFromFasta(arrayOfFasta[i]);
    title = getTitleFromFasta(arrayOfFasta[i]);
    newDna = removeNonDna(newDna);

    //create restriction fragment objects from sequences and add them to array
    restrictionFragment = new RestrictionFragment(
      theDocument.forms[0].elements[4].options[
        theDocument.forms[0].elements[4].selectedIndex
      ].value,
      title,
      newDna,
      1,
      newDna.length,
      "sequence start",
      "sequence end",
      newDna.length
    );
    restrictionFragments.push(restrictionFragment);
  }

  restrictionFragments = digest(
    restrictionFragments,
    theDocument.forms[0].elements[5].options[
      theDocument.forms[0].elements[5].selectedIndex
    ].value
  );
  restrictionFragments = digest(
    restrictionFragments,
    theDocument.forms[0].elements[6].options[
      theDocument.forms[0].elements[6].selectedIndex
    ].value
  );
  restrictionFragments = digest(
    restrictionFragments,
    theDocument.forms[0].elements[7].options[
      theDocument.forms[0].elements[7].selectedIndex
    ].value
  );

  restrictionFragments.sort(restrictionFragmentSorter);

  openWindow("Restriction Digest");
  openPre();
  for (var i = 0; i < restrictionFragments.length; i++) {
    restrictionFragments[i].correctPositions();
    restrictionFragments[i].writeFragment(
      theDocument.forms[0].elements[4].options[
        theDocument.forms[0].elements[4].selectedIndex
      ].value
    );
  }
  closePre();
  closeWindow();
  return true;
}

function digest(arrayOfRestrictionFragments, enzyme) {
  if (enzyme == "") {
    return arrayOfRestrictionFragments;
  }

  var newFragments = new Array();
  var positions = new Array();
  var matchPosition = 0;
  var matchExp = enzyme.match(/\/.+\//).toString();
  matchExp = matchExp.replace(/\//g, "");
  matchExp = new RegExp(matchExp, "gi");
  var cutDistance = parseInt(
    enzyme
      .match(/\)\D*\d+/)
      .toString()
      .replace(/\)\D*/, "")
  );
  var enzymeName = enzyme
    .match(/\([^\(]+\)/)
    .toString()
    .replace(/\(|\)/g, "");
  var matchArray;
  var wasCut = false;
  var restrictionFragmentOne;
  var restrictionFragmentTwo;

  var previousCutPosition;
  var previousEnzyme;
  var previousStartPosition;

  var startRestrictionFragment;

  var lookAhead = 50;
  var lowerLimit;
  var upperLimit;
  var shiftValue;

  //go through each restriction fragment.
  for (var i = 0; i < arrayOfRestrictionFragments.length; i++) {
    if (arrayOfRestrictionFragments[i].topology == "circular") {
      shiftValue = arrayOfRestrictionFragments[i].sequence.substring(
        0,
        lookAhead
      ).length;
      var extendedSequence =
        arrayOfRestrictionFragments[i].sequence.substring(
          arrayOfRestrictionFragments[i].sequence.length - lookAhead,
          arrayOfRestrictionFragments[i].sequence.length
        ) +
        arrayOfRestrictionFragments[i].sequence +
        arrayOfRestrictionFragments[i].sequence.substring(0, lookAhead);
      lowerLimit = 0 + shiftValue;
      upperLimit = arrayOfRestrictionFragments[i].sequence.length + shiftValue;

      while ((matchArray = matchExp.exec(extendedSequence))) {
        matchPosition = matchExp.lastIndex;
        matchPosition = matchPosition - cutDistance;
        if (matchPosition >= lowerLimit && matchPosition < upperLimit) {
          positions.push(matchPosition - shiftValue);
          wasCut = true;
        }
        matchExp.lastIndex = matchExp.lastIndex - RegExp.lastMatch.length + 1;
      }
    } else {
      while (
        (matchArray = matchExp.exec(arrayOfRestrictionFragments[i].sequence))
      ) {
        matchPosition = matchExp.lastIndex;
        matchPosition = matchPosition - cutDistance;
        positions.push(matchPosition);
        wasCut = true;
        matchExp.lastIndex = matchExp.lastIndex - RegExp.lastMatch.length + 1;
      }
    }

    if (!wasCut) {
      newFragments.push(arrayOfRestrictionFragments[i]);
    } else {
      wasCut = false;

      //now go through positions
      previousCutPosition = 0;
      previousEnzyme = arrayOfRestrictionFragments[i].startEnzyme;
      previousStartPosition = arrayOfRestrictionFragments[i].start;
      for (var j = 0; j < positions.length; j++) {
        if (arrayOfRestrictionFragments[i].topology == "circular") {
          arrayOfRestrictionFragments[i].topology = "linear";
          startRestrictionFragment = new RestrictionFragment(
            "linear",
            arrayOfRestrictionFragments[i].sourceName,
            arrayOfRestrictionFragments[i].sequence.substring(
              previousCutPosition,
              positions[j]
            ),
            previousStartPosition,
            previousStartPosition +
              arrayOfRestrictionFragments[i].sequence.substring(
                previousCutPosition,
                positions[j]
              ).length -
              1,
            previousEnzyme,
            enzymeName,
            arrayOfRestrictionFragments[i].originalLength
          );
        } else {
          restrictionFragmentOne = new RestrictionFragment(
            "linear",
            arrayOfRestrictionFragments[i].sourceName,
            arrayOfRestrictionFragments[i].sequence.substring(
              previousCutPosition,
              positions[j]
            ),
            previousStartPosition,
            previousStartPosition +
              arrayOfRestrictionFragments[i].sequence.substring(
                previousCutPosition,
                positions[j]
              ).length -
              1,
            previousEnzyme,
            enzymeName,
            arrayOfRestrictionFragments[i].originalLength
          );
          newFragments.push(restrictionFragmentOne);
        }

        if (j == positions.length - 1) {
          if (startRestrictionFragment == null) {
            restrictionFragmentTwo = new RestrictionFragment(
              "linear",
              arrayOfRestrictionFragments[i].sourceName,
              arrayOfRestrictionFragments[i].sequence.substring(
                positions[j],
                arrayOfRestrictionFragments[i].sequence.length
              ),
              previousStartPosition +
                arrayOfRestrictionFragments[i].sequence.substring(
                  previousCutPosition,
                  positions[j]
                ).length,
              arrayOfRestrictionFragments[i].stop,
              enzymeName,
              arrayOfRestrictionFragments[i].stopEnzyme,
              arrayOfRestrictionFragments[i].originalLength
            );
            newFragments.push(restrictionFragmentTwo);
          } else {
            restrictionFragmentTwo = new RestrictionFragment(
              "linear",
              arrayOfRestrictionFragments[i].sourceName,
              arrayOfRestrictionFragments[i].sequence.substring(
                positions[j],
                arrayOfRestrictionFragments[i].sequence.length
              ) + startRestrictionFragment.sequence,
              previousStartPosition +
                arrayOfRestrictionFragments[i].sequence.substring(
                  previousCutPosition,
                  positions[j]
                ).length,
              startRestrictionFragment.stop,
              enzymeName,
              startRestrictionFragment.stopEnzyme,
              arrayOfRestrictionFragments[i].originalLength
            );

            newFragments.push(restrictionFragmentTwo);
          }
        }

        previousCutPosition = positions[j];
        previousEnzyme = enzymeName;
        previousStartPosition =
          arrayOfRestrictionFragments[i].start + positions[j];
      }
    }
    positions = new Array();
    startRestrictionFragment = null;
  }
  return newFragments;
}

function restrictionFragmentSorter(a, b) {
  if (a.sequence.length < b.sequence.length) {
    return 1;
  }
  if (a.sequence.length > b.sequence.length) {
    return -1;
  }
  if (a.sequence.length == b.sequence.length) {
    if (a.start < b.start) {
      return -1;
    }
    if (a.start > b.start) {
      return 1;
    } else {
      // a must be equal to b
      return 0;
    }
  }
}

//RestrictionFragment writeFragment() method
function writeFragment(topology) {
  if (this.topology == "linear") {
    outputWindow.document.write(
      "&gt;" +
        this.sequence.length +
        " bp linear fragment from " +
        topology +
        " parent " +
        this.sourceName +
        ", base " +
        this.start +
        " to base " +
        this.stop +
        " (" +
        this.startEnzyme +
        " - " +
        this.stopEnzyme +
        ").\n"
    );
  } else {
    outputWindow.document.write(
      "&gt;" +
        this.sequence.length +
        " bp circular molecule from " +
        topology +
        " parent " +
        this.sourceName +
        ".\n"
    );
  }
  outputWindow.document.write(addReturns(this.sequence) + "\n\n");
}

//RestrictionFragment correctPositions() method
function correctPositions() {
  if (this.start > this.originalLength) {
    this.start = this.start - this.originalLength;
  }
  if (this.stop > this.originalLength) {
    this.stop = this.stop - this.originalLength;
  }
  if (this.stop == 0) {
    this.stop = this.originalLength;
  }
}

//RestrictionFragment class
function RestrictionFragment(
  topology,
  sourceName,
  sequence,
  start,
  stop,
  startEnzyme,
  stopEnzyme,
  originalLength
) {
  this.topology = topology;
  this.sourceName = sourceName;
  this.sequence = sequence;
  this.start = start;
  this.stop = stop;
  this.startEnzyme = startEnzyme;
  this.stopEnzyme = stopEnzyme;
  this.originalLength = originalLength;
}

//create and throw away a prototype object
new RestrictionFragment("", "", "", 0, 0, "", "", 0);

// define object methods
RestrictionFragment.prototype.writeFragment = writeFragment;
RestrictionFragment.prototype.correctPositions = correctPositions;
