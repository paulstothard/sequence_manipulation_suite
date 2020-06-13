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

function threeToOne(theDocument) {
  var newProtein = "";
  var maxInput = 100000000;

  if (testScript() == false) {
    return false;
  }

  if (
    checkFormElement(theDocument.forms[0].elements[0]) == false ||
    checkTextLength(theDocument.forms[0].elements[0].value, maxInput) == false
  ) {
    return false;
  }

  openWindow("Three to One");
  openPre();
  var arrayOfFasta = getArrayOfFasta(theDocument.forms[0].elements[0].value);

  for (var i = 0; i < arrayOfFasta.length; i++) {
    newProtein = getTripletSequenceFromFasta(arrayOfFasta[i]);
    title = getFastaTitleFromTriplets(arrayOfFasta[i]);
    newProtein = filterTriplets(newProtein);

    outputWindow.document.write(
      getInfoFromTitleAndSequenceTriplets(title, newProtein)
    );

    writeThreeToOne(newProtein);

    outputWindow.document.write("\n\n");
  }

  closePre();
  closeWindow();
  return true;
}

function writeThreeToOne(proteinSequence) {
  proteinSequence = proteinSequence.replace(/(.)(.)(.)/g, function (
    str,
    p1,
    p2,
    p3,
    offset,
    s
  ) {
    return p1.toUpperCase() + p2.toLowerCase() + p3.toLowerCase();
  });
  proteinSequence = proteinSequence.replace(/Ala/g, " A ");
  proteinSequence = proteinSequence.replace(/Asx/g, " B ");
  proteinSequence = proteinSequence.replace(/Cys/g, " C ");
  proteinSequence = proteinSequence.replace(/Asp/g, " D ");
  proteinSequence = proteinSequence.replace(/Glu/g, " E ");
  proteinSequence = proteinSequence.replace(/Phe/g, " F ");
  proteinSequence = proteinSequence.replace(/Gly/g, " G ");
  proteinSequence = proteinSequence.replace(/His/g, " H ");
  proteinSequence = proteinSequence.replace(/Ile/g, " I ");
  proteinSequence = proteinSequence.replace(/Lys/g, " K ");
  proteinSequence = proteinSequence.replace(/Leu/g, " L ");
  proteinSequence = proteinSequence.replace(/Met/g, " M ");
  proteinSequence = proteinSequence.replace(/Asn/g, " N ");
  proteinSequence = proteinSequence.replace(/Pro/g, " P ");
  proteinSequence = proteinSequence.replace(/Gln/g, " Q ");
  proteinSequence = proteinSequence.replace(/Arg/g, " R ");
  proteinSequence = proteinSequence.replace(/Ser/g, " S ");
  proteinSequence = proteinSequence.replace(/Thr/g, " T ");
  proteinSequence = proteinSequence.replace(/Val/g, " V ");
  proteinSequence = proteinSequence.replace(/Trp/g, " W ");
  proteinSequence = proteinSequence.replace(/Xaa/g, " X ");
  proteinSequence = proteinSequence.replace(/Tyr/g, " Y ");
  proteinSequence = proteinSequence.replace(/Glx/g, " Z ");
  proteinSequence = proteinSequence.replace(/\*\*\*/g, " * ");

  proteinSequence = proteinSequence.replace(/\s/g, "");

  outputWindow.document.write(addReturns(proteinSequence));
  return true;
}

function filterTriplets(tripletSequence) {
  tripletSequence = tripletSequence.replace(/\s|\d/gi, "");
  return tripletSequence;
}

function getFastaTitleFromTriplets(tripletSequence) {
  fastaSequenceTitle = "Untitled";
  if (tripletSequence.search(/\>[^\f\n\r]+[\f\n\r]/) != -1) {
    fastaSequenceTitle = tripletSequence
      .match(/\>[^\f\n\r]+[\f\n\r]/, "")
      .toString();
    fastaSequenceTitle = fastaSequenceTitle.replace(/\>|[\f\n\r]/g, "");
    fastaSequenceTitle = filterFastaTitle(fastaSequenceTitle);
  }
  return fastaSequenceTitle;
}

function getTripletSequenceFromFasta(tripletSequence) {
  if (tripletSequence.search(/\>[^\f\n\r]+[\f\n\r]/) != -1) {
    tripletSequence = tripletSequence.replace(/\>[^\f\n\r]+[\f\n\r]/, "");
  }
  return tripletSequence;
}

function getInfoFromTitleAndSequenceTriplets(fastaSequenceTitle, sequence) {
  var stringToReturn = "&gt;results for sequence ";
  if (fastaSequenceTitle.search(/[^\s]/) != -1) {
    stringToReturn = stringToReturn + '"' + fastaSequenceTitle + '"';
  }
  stringToReturn =
    stringToReturn + ' starting "' + sequence.substring(0, 12) + '"';
  return stringToReturn + "\n";
}
