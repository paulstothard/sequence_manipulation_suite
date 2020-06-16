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

function proteinIep(theDocument) {
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

  openWindow("Protein Isoelectric Point");
  var arrayOfFasta = getArrayOfFasta(theDocument.forms[0].elements[0].value);

  for (var i = 0; i < arrayOfFasta.length; i++) {
    newProtein = getSequenceFromFasta(arrayOfFasta[i]);
    title = getTitleFromFasta(arrayOfFasta[i]);
    newProtein = removeNonProteinStrict(newProtein);

    outputWindow.document.write(getInfoFromTitleAndSequence(title, newProtein));

    writeProtIep(
      newProtein,
      theDocument.forms[0].elements[4].options[
        theDocument.forms[0].elements[4].selectedIndex
      ].value,
      theDocument.forms[0].elements[5].options[
        theDocument.forms[0].elements[5].selectedIndex
      ].value,
      theDocument.forms[0].elements[6].options[
        theDocument.forms[0].elements[6].selectedIndex
      ].value
    );

    outputWindow.document.write("<br />\n<br />\n");
  }
  closeWindow();
  return true;
}

function writeProtIep(proteinSequence, copies, fusion, pKSet) {
  //calculates pI of protein.
  var pH = 7.0;
  var step = 3.5;
  var charge = 0.0;
  var last_charge = 0.0;

  copies = parseInt(copies);
  for (var j = 0; j < copies; j++) {
    proteinSequence = proteinSequence + fusion;
  }

  var N_term_pK;
  var K_pK;
  var R_pK;
  var H_pK;
  var D_pK;
  var E_pK;
  var C_pK;
  var Y_pK;
  var C_term_pK;

  if (pKSet.toLowerCase() == "dtaselect") {
    //pK values from DTASelect
    N_term_pK = 8.0;
    K_pK = 10.0;
    R_pK = 12.0;
    H_pK = 6.5;
    D_pK = 4.4;
    E_pK = 4.4;
    C_pK = 8.5;
    Y_pK = 10.0;
    C_term_pK = 3.1;
  } else {
    //pK values from EMBOSS
    N_term_pK = 8.6;
    K_pK = 10.8;
    R_pK = 12.5;
    H_pK = 6.5;
    D_pK = 3.9;
    E_pK = 4.1;
    C_pK = 8.5;
    Y_pK = 10.1;
    C_term_pK = 3.6;
  }

  var K_count = 0;
  if (proteinSequence.search(/k/i) != -1) {
    K_count = proteinSequence.match(/k/gi).length;
  }

  var R_count = 0;
  if (proteinSequence.search(/r/i) != -1) {
    R_count = proteinSequence.match(/r/gi).length;
  }

  var H_count = 0;
  if (proteinSequence.search(/h/i) != -1) {
    H_count = proteinSequence.match(/h/gi).length;
  }

  var D_count = 0;
  if (proteinSequence.search(/d/i) != -1) {
    D_count = proteinSequence.match(/d/gi).length;
  }

  var E_count = 0;
  if (proteinSequence.search(/e/i) != -1) {
    E_count = proteinSequence.match(/e/gi).length;
  }

  var C_count = 0;
  if (proteinSequence.search(/c/i) != -1) {
    C_count = proteinSequence.match(/c/gi).length;
  }

  var Y_count = 0;
  if (proteinSequence.search(/y/i) != -1) {
    Y_count = proteinSequence.match(/y/gi).length;
  }

  while (1) {
    charge =
      partial_charge(N_term_pK, pH) +
      K_count * partial_charge(K_pK, pH) +
      R_count * partial_charge(R_pK, pH) +
      H_count * partial_charge(H_pK, pH) -
      D_count * partial_charge(pH, D_pK) -
      E_count * partial_charge(pH, E_pK) -
      C_count * partial_charge(pH, C_pK) -
      Y_count * partial_charge(pH, Y_pK) -
      partial_charge(pH, C_term_pK);

    if (charge.toFixed(2) == (last_charge * 100).toFixed(2)) {
      break;
    }

    if (charge > 0) {
      pH = pH + step;
    } else {
      pH = pH - step;
    }

    step = step / 2;

    last_charge = charge;
  }

  pH = pH.toFixed(2);
  outputWindow.document.write("pH " + pH);

  return true;
}

function partial_charge(first, second) {
  var charge = Math.pow(10, first - second);
  return charge / (charge + 1);
}
