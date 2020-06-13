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

function emblTrans(theDocument) {
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

  openWindow("EMBL Trans Extractor");
  openPre();
  emblTransExtract(theDocument.forms[0].elements[0].value);
  closePre();
  closeWindow();
  return true;
}
function emblTransExtract(emblFile) {
  var title;
  emblFile = "_" + emblFile + "_";
  var recordArray = emblFile.split(/ID\s\s\s[^\f\n\r]*/);
  for (var i = 1; i < recordArray.length; i++) {
    var mainArray = emblFile.split(
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
    var featureArray = mainArray[1].split(/[\f\n\r]FT {3,12}\b/);
    showFeatureTrans(featureArray);
  }
  return true;
}

function showFeatureTrans(arrayOfFeatures) {
  var featureTitle = "";
  var theTitle = "";
  var removedTitle = "";
  var firstQualifier = "";
  var translation = "";
  var translationFound = false;
  for (var i = 1; i < arrayOfFeatures.length; i++) {
    if (arrayOfFeatures[i].search(/\/translation/) != -1) {
      arrayOfFeatures[i] = arrayOfFeatures[i].replace(/[\[\]\*]/g, "");
      featureTitle = arrayOfFeatures[i].match(/[^ \f\n\r\t\v]+ /).toString();
      theTitle = new RegExp(featureTitle);
      removedTitle = arrayOfFeatures[i].replace(theTitle, "");
      firstQualifier = arrayOfFeatures[i].match(/\/[^\f\n\r]+/).toString();
      outputWindow.document.write(
        "&gt;" +
          filterFastaTitle(featureTitle) +
          filterFastaTitle(firstQualifier) +
          "\n"
      );
      translation = arrayOfFeatures[i].match(/\/translation="[^"]+"/);
      translation = translation.toString();
      translation = translation.replace(/\/translation/, "");
      translation = translation.replace(/^FT\s+/gm, "");
      translation = removeNonProtein(translation);
      translation = addReturns(translation);
      outputWindow.document.write(translation);
      translationFound = true;
      outputWindow.document.write("\n\n");
    }
  }
  if (translationFound == false) {
    outputWindow.document.write("No translations were found.\n");
  }
  return true;
}
