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

function dnaMw(theDocument) {
  var newDna = "";
  var maxInput = 200000000;

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

  openWindow("DNA Molecular Weight");
  var arrayOfFasta = getArrayOfFasta(theDocument.forms[0].elements[0].value);

  for (var i = 0; i < arrayOfFasta.length; i++) {
    newDna = getSequenceFromFasta(arrayOfFasta[i]);
    title = getTitleFromFasta(arrayOfFasta[i]);
    newDna = _removeNonPrimer(newDna);

    outputWindow.document.write(getInfoFromTitleAndSequence(title, newDna));

    writeDnaMw(
      newDna,
      theDocument.forms[0].elements[4].options[
        theDocument.forms[0].elements[4].selectedIndex
      ].value,
      theDocument.forms[0].elements[5].options[
        theDocument.forms[0].elements[5].selectedIndex
      ].value
    );

    outputWindow.document.write("<br />\n<br />\n");
  }
  closeWindow();
  return true;
}

function writeDnaMw(sequence, strandType, topology) {
  //calculates molecular weight of DNA.
  //ligation removes OH
  var OH = 17.01;
  var result = 0;

  if (strandType == "single") {
    var mw_direct_strand = _molecularWeight(sequence);
    if (mw_direct_strand.length == 1) {
      var mw = parseFloat(mw_direct_strand[0]);
      if (topology == "circular") {
        mw = mw - OH;
      }
      mw = mw.toFixed(2);
      outputWindow.document.write(mw + " Da");
    } else if (mw_direct_strand.length == 2) {
      var mw_lower = parseFloat(mw_direct_strand[0]);
      var mw_upper = parseFloat(mw_direct_strand[1]);
      if (topology == "circular") {
        mw_lower = mw_lower - OH;
        mw_upper = mw_upper - OH;
      }
      mw_lower = mw_lower.toFixed(2);
      mw_upper = mw_upper.toFixed(2);
      outputWindow.document.write(mw_lower + " to " + mw_upper + " Da");
    }
  } else if (strandType == "double") {
    var mw_direct_strand = _molecularWeight(sequence);
    var mw_reverse_strand = _molecularWeight(reverse(complement(sequence)));
    if (mw_direct_strand.length == 1 && mw_reverse_strand.length == 1) {
      var mw_direct = parseFloat(mw_direct_strand[0]);
      var mw_reverse = parseFloat(mw_reverse_strand[0]);
      if (topology == "circular") {
        mw_direct = mw_direct - OH;
        mw_reverse = mw_reverse - OH;
      }
      var mw = mw_direct + mw_reverse;
      mw = mw.toFixed(2);
      outputWindow.document.write(mw + " Da");
    } else if (mw_direct_strand.length == 2 && mw_reverse_strand.length == 2) {
      var mw_direct_lower = parseFloat(mw_direct_strand[0]);
      var mw_reverse_lower = parseFloat(mw_reverse_strand[0]);
      var mw_direct_upper = parseFloat(mw_direct_strand[1]);
      var mw_reverse_upper = parseFloat(mw_reverse_strand[1]);
      if (topology == "circular") {
        mw_direct_lower = mw_direct_lower - OH;
        mw_reverse_lower = mw_reverse_lower - OH;

        mw_direct_upper = mw_direct_upper - OH;
        mw_reverse_upper = mw_reverse_upper - OH;
      }
      var mw_lower = mw_direct_lower + mw_reverse_lower;
      var mw_upper = mw_direct_upper + mw_reverse_upper;
      mw_lower = mw_lower.toFixed(2);
      mw_upper = mw_upper.toFixed(2);
      outputWindow.document.write(mw_lower + " to " + mw_upper + " Da");
    }
  }

  return true;
}

function _containsOnlyNonDegenerates(sequence) {
  if (sequence.search(/[^gatc]/i) == -1) {
    return true;
  }
  return false;
}

function _molecularWeight(sequence) {
  if (_containsOnlyNonDegenerates(sequence)) {
    return _molecularWeightNonDegen(sequence);
  } else {
    return _molecularWeightDegen(sequence);
  }
}

function _molecularWeightNonDegen(sequence) {
  var results = new Array();
  results[0] = _mw(sequence);
  return results;
}

function _mw(sequence) {
  //DNA molecular weight for linear strand of DNA with a 5' monophosphate
  var g = _getBaseCount(sequence, "g");
  var a = _getBaseCount(sequence, "a");
  var t = _getBaseCount(sequence, "t");
  var c = _getBaseCount(sequence, "c");
  return g * 329.21 + a * 313.21 + t * 304.2 + c * 289.18 + 17.01;
}

function _molecularWeightDegen(sequence) {
  var lowerBoundsSequence = sequence;
  var upperBoundsSequence = sequence;

  //replace all other degenerates with lightest base possible in lowerBoundsSequence
  lowerBoundsSequence = lowerBoundsSequence.replace(/r/gi, "a");
  lowerBoundsSequence = lowerBoundsSequence.replace(/y/gi, "c");
  lowerBoundsSequence = lowerBoundsSequence.replace(/s/gi, "c");
  lowerBoundsSequence = lowerBoundsSequence.replace(/w/gi, "t");
  lowerBoundsSequence = lowerBoundsSequence.replace(/k/gi, "t");
  lowerBoundsSequence = lowerBoundsSequence.replace(/m/gi, "c");
  lowerBoundsSequence = lowerBoundsSequence.replace(/b/gi, "c");
  lowerBoundsSequence = lowerBoundsSequence.replace(/d/gi, "t");
  lowerBoundsSequence = lowerBoundsSequence.replace(/h/gi, "c");
  lowerBoundsSequence = lowerBoundsSequence.replace(/v/gi, "c");
  lowerBoundsSequence = lowerBoundsSequence.replace(/n/gi, "c");

  //replace all other degenerates with heaviest base possible in upperBoundsSequence
  upperBoundsSequence = upperBoundsSequence.replace(/r/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/y/gi, "t");
  upperBoundsSequence = upperBoundsSequence.replace(/s/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/w/gi, "a");
  upperBoundsSequence = upperBoundsSequence.replace(/k/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/m/gi, "a");
  upperBoundsSequence = upperBoundsSequence.replace(/b/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/d/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/h/gi, "a");
  upperBoundsSequence = upperBoundsSequence.replace(/v/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/n/gi, "g");

  var results = new Array();
  results[0] = _molecularWeightNonDegen(lowerBoundsSequence);
  results[1] = _molecularWeightNonDegen(upperBoundsSequence);
  return results;
}

function _getBaseCount(sequence, base) {
  var basePattern = new RegExp(base, "gi");
  if (sequence.search(basePattern) != -1) {
    return sequence.match(basePattern).length;
  } else {
    return 0;
  }
}

function _removeNonPrimer(sequence) {
  sequence.replace(/u/g, "t");
  sequence.replace(/U/g, "T");
  return sequence.replace(/[^gatcryswkmbdhvnGATCRYSWKMBDHVN]/g, "");
}
