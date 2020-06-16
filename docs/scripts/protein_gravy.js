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

function proteinGravy(theDocument) {
  var newProtein = "";
  var title = "";
  var maxInput = 500000000;

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

  openWindow("Protein GRAVY");
  var arrayOfFasta = getArrayOfFasta(theDocument.forms[0].elements[0].value);

  for (var i = 0; i < arrayOfFasta.length; i++) {
    newProtein = getSequenceFromFasta(arrayOfFasta[i]);

    title = getTitleFromFasta(arrayOfFasta[i]);

    newProtein = removeNonProtein(newProtein);

    outputWindow.document.write(getInfoFromTitleAndSequence(title, newProtein));

    outputWindow.document.write(getProteinGravy(newProtein));

    outputWindow.document.write("<br />\n<br />\n");
  }

  closeWindow();
  return true;
}

function getProteinGravy(sequence) {
  sequence = sequence.toLowerCase();
  var gravyResult = 0;
  //The GRAVY value for a peptide or protein is calculated as the sum of hydropathy values [9]
  //of all the amino acids, divided by the number of residues in the sequence.
  var gravyValues = _getGravyHash();
  for (var i = 0; i < sequence.length; i++) {
    gravyResult = gravyResult + gravyValues[sequence.charAt(i)];
  }
  if (sequence.length > 0) {
    gravyResult = gravyResult / sequence.length;
  } else {
    return "The sequence is too short";
  }
  return gravyResult.toFixed(3);
}

function _getGravyHash() {
  //Author(s): Kyte J., Doolittle R.F.
  //Reference: J. Mol. Biol. 157:105-132(1982).
  var hash = {};
  hash["a"] = 1.8;
  hash["r"] = -4.5;
  hash["n"] = -3.5;
  hash["d"] = -3.5;
  hash["c"] = 2.5;
  hash["q"] = -3.5;
  hash["e"] = -3.5;
  hash["g"] = -0.4;
  hash["h"] = -3.2;
  hash["i"] = 4.5;
  hash["l"] = 3.8;
  hash["k"] = -3.9;
  hash["m"] = 1.9;
  hash["f"] = 2.8;
  hash["p"] = -1.6;
  hash["s"] = -0.8;
  hash["t"] = -0.7;
  hash["w"] = -0.9;
  hash["y"] = -1.3;
  hash["v"] = 4.2;
  return hash;
}
