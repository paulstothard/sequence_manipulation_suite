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

function proteinMw(theDocument) {
  var newProtein = "";
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

  //the weights below have water subtracted
  var arrayOfMw = [
    "/A/ (A)71.08",
    "/C/ (C)103.14",
    "/D/ (D)115.09",
    "/E/ (E)129.12",
    " /F/ (F)147.18",
    "/G/ (G)57.06",
    "/H/ (H)137.15",
    "/I/ (I)113.17",
    "/K/ (K)128.18",
    "/L/ (L)113.17",
    "/M/ (M)131.21",
    "/N/ (N)114.11",
    "/P/ (P)97.12",
    "/Q/ (Q)128.41",
    "/R/ (R)156.20",
    "/S/ (S)87.08",
    "/T/ (T)101.11",
    "/V/ (V)99.14",
    "/W/ (W)186.21",
    "/Y/ (Y)163.18",
  ];

  openWindow("Protein Molecular Weight");
  var arrayOfFasta = getArrayOfFasta(theDocument.forms[0].elements[0].value);

  for (var i = 0; i < arrayOfFasta.length; i++) {
    newProtein = getSequenceFromFasta(arrayOfFasta[i]);
    title = getTitleFromFasta(arrayOfFasta[i]);
    newProtein = removeNonProteinStrict(newProtein);

    outputWindow.document.write(getInfoFromTitleAndSequence(title, newProtein));

    writeProtMw(
      newProtein,
      arrayOfMw,
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

function writeProtMw(proteinSequence, arrayOfMw, copies, fusion) {
  //calculates molecular weight of protein.
  var water = 18.015;
  var result = 0;
  copies = parseInt(copies);
  for (var j = 0; j < copies; j++) {
    proteinSequence = proteinSequence + fusion;
  }
  for (var j = 0; j < arrayOfMw.length; j++) {
    var tempNumber = 0;
    var matchExp = arrayOfMw[j].match(/\/[^\/]+\//) + "gi";
    matchExp = eval(matchExp);
    if (proteinSequence.search(matchExp) != -1) {
      tempNumber = proteinSequence.match(matchExp).length;
    }
    result =
      result +
      tempNumber * parseFloat(arrayOfMw[j].match(/[\d\.]+/).toString());
  }

  if (result == 0) {
    outputWindow.document.write(result + " kDa");
  } else {
    result = result + water; //add the weight of water for the ends of the protein.
    result = result / 1000; //convert to kilodaltons.
    result = result.toFixed(2);
    outputWindow.document.write(result + " kDa");
  }
  return true;
}
