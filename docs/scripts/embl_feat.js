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

function emblFeat(theDocument) {
  var maxInput = 200000000;

  if (testScript() == false) {
    return false;
  }

  if (
    checkFormElement(theDocument.forms[0].elements[0]) == false ||
    verifyEmblFeat(theDocument.forms[0].elements[0].value) == false ||
    checkTextLength(theDocument.forms[0].elements[0].value, maxInput) == false
  ) {
    return false;
  }

  openWindow("EMBL Feature Extractor");
  openPre();
  emblFeatExtract(
    theDocument.forms[0].elements[0].value,
    theDocument.forms[0].elements[4].options[
      theDocument.forms[0].elements[4].selectedIndex
    ].value
  );
  closePre();
  closeWindow();
  return true;
}

function emblFeatExtract(emblFile, outputType) {
  var title;
  emblFile = "_" + emblFile + "_";
  var recordArray = emblFile.split(/ID\s\s\s[^\f\n\r]*/);
  for (var i = 1; i < recordArray.length; i++) {
    var mainArray = recordArray[i].split(
      /[\f\n\r]\s*FH   Key             Location\/Qualifiers[\f\n\r]+\s*FH|[\f\n\r]\s*XX[\s]*[\f\n\r]\s*SQ[^\f\n\r]*/
    );
    if (mainArray[0].search(/[\f\n\r]\s*DE[^\f\n\r]+/) != -1) {
      title = mainArray[0]
        .match(/[\f\n\r]\s*DE[^\f\n\r]+/)
        .toString()
        .replace(/[\f\n\r]\s*DE\s*/, "");
    } else {
      title = "Untitled";
    }
    title = filterFastaTitle(title.replace(/[\f\n\r\t]+$/g, "")) + "\n";
    var dnaSequenceArray = mainArray[2].split(/\/{2}/);
    outputWindow.document.write(title + "\n");
    if (dnaSequenceArray.length == 1) {
      alert("The entire EMBL file may not have been processed.");
      outputWindow.focus();
    }
    var dnaSequence = removeNonDna(dnaSequenceArray[0]);
    var featureArray = mainArray[1].split(/[\f\n\r]FT {3,12}\b/);
    prepareFeatures(featureArray, dnaSequence, outputType);
  }
  return true;
}

function prepareFeatures(arrayOfFeatures, dnaSequence, outputType) {
  var featureTitle = "";
  var theTitle = "";
  var removedTitle = "";
  var firstQualifier = "";
  var position = "";
  var positionNoSpace = "";
  var featureFound = false;
  for (var i = 1; i < arrayOfFeatures.length; i++) {
    arrayOfFeatures[i] = arrayOfFeatures[i].replace(/[\[\]\*]/g, "");
    featureTitle = arrayOfFeatures[i].match(/[^ \f\n\r\t\v]+ /).toString();
    theTitle = new RegExp(featureTitle);
    removedTitle = arrayOfFeatures[i].replace(theTitle, "");
    if (arrayOfFeatures[i].search(/\/[^\f\n\r]+/) != -1) {
      firstQualifier = arrayOfFeatures[i].match(/\/[^\f\n\r]+/).toString();
    } else {
      firstQualifier = "/no qualifier supplied";
    }
    position = removedTitle.split(/\//);
    positionNoSpace = position[0].replace(/\s{2,}/g, " ");
    outputWindow.document.write(
      "&gt;" +
        filterFastaTitle(featureTitle) +
        filterFastaTitle(firstQualifier) +
        "\n"
    );
    printFeature(positionNoSpace, dnaSequence, outputType);
    featureFound = true;
    outputWindow.document.write("\n\n");
  }
  if (featureFound == false) {
    outputWindow.document.write(
      "There were no features found or there was a problem reading the feature information."
    );
  }
  return true;
}

function printFeature(featurePos, dnaSequence, outputType) {
  var feature;
  featurePos = featurePos.replace(/<|>/g, "");
  featurePos = featurePos.replace(/FT/gi, "");
  if (
    featurePos.search(/[^a-z\d\.\(\)\,\s]/) != -1 ||
    featurePos.search(/one/) != -1 ||
    featurePos.search(/order/) != -1
  ) {
    outputWindow.document.write(
      "This feature specifies a sequence that cannot be represented:\n"
    );
    outputWindow.document.write(featurePos);
  } else {
    var newFeaturePos = featurePos.replace(/\)/g, "");
    if (newFeaturePos.search(/complement/) != -1) {
      feature = new Feature("complement");
    } else {
      feature = new Feature("direct");
    }

    var pairString = newFeaturePos;
    var pairArray = pairString.split(/\,/);

    if (newFeaturePos.search(/complement/) != -1) {
      pairArray.reverse();
    }

    var digitArray = new Array();
    var realStart = 0;
    var realStop = 0;
    var dnaToAdd = "";

    for (var j = 0; j < pairArray.length; j++) {
      digitArray = pairArray[j].split(/\.\./);
      if (digitArray.length == 1) {
        digitArray[1] = digitArray[0];
      }
      realStart = digitArray[0];
      realStop = digitArray[1];
      realStop = realStop.replace(/\D/g, "");
      realStart = realStart.replace(/\D/g, "");
      if (realStart.search(/\d/) == -1 || realStop.search(/\d/) == -1) {
        outputWindow.document.write(
          "There was a problem with this feature (one of the range values was missing)."
        );
        return true;
      }
      realStart = parseInt(realStart) - 1;
      realStop = parseInt(realStop);
      if (realStart > realStop) {
        outputWindow.document.write(
          "There was a problem with this feature (the end position was before the start position)."
        );
        return true;
      }
      if (realStart > dnaSequence.length || realStop > dnaSequence.length) {
        outputWindow.document.write(
          "The entire EMBL file was not processed, so this feature cannot be properly shown."
        ); //to handle case where web browser limits the number of characters that can be entered as input.
        return true;
      } else {
        if (outputType == "separated") {
          feature.addFragment(dnaSequence.substring(realStart, realStop));
        } else {
          feature.addFragment(
            dnaSequence.substring(feature.lastAdded, realStart)
          );
          feature.addFragment(
            dnaSequence.substring(realStart, realStop).toUpperCase()
          );
          feature.lastAdded = realStop;
        }
      }
    }
    feature.writeFeature();
  }
  return true;
}

//class Feature method writeFeature()
function writeFeature() {
  if (this.strand == "complement") {
    outputWindow.document.write(
      addReturns(reverse(complement(this.fragments.join(""))))
    );
  } else {
    outputWindow.document.write(addReturns(this.fragments.join("")));
  }
}

//class Feature method addFragment()
function addFragment(sequence) {
  this.fragments.push(sequence);
}

//class Feature
function Feature(strand) {
  this.strand = strand;
  this.fragments = new Array();
  this.lastAdded = 0;
}

//create and throw away a prototype object
new Feature();

// define object methods
Feature.prototype.writeFeature = writeFeature;
Feature.prototype.addFragment = addFragment;
