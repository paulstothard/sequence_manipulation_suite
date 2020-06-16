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

function oneToThree(theDocument) {
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

  openWindow("One to Three");
  openPre();
  var arrayOfFasta = getArrayOfFasta(theDocument.forms[0].elements[0].value);

  for (var i = 0; i < arrayOfFasta.length; i++) {
    newProtein = getSequenceFromFasta(arrayOfFasta[i]);
    title = getTitleFromFasta(arrayOfFasta[i]);
    newProtein = removeNonProteinAllowDegen(newProtein);

    outputWindow.document.write(
      getFastaTitleFromTitleAndSequence(title, newProtein)
    );

    writeOneToThree(newProtein);

    outputWindow.document.write("\n\n");
  }

  closePre();
  closeWindow();
  return true;
}

function writeOneToThree(proteinSequence) {
  proteinSequence = proteinSequence.toLowerCase();
  proteinSequence = proteinSequence.replace(/(.)/g, function (
    str,
    p1,
    offset,
    s
  ) {
    return " " + p1 + " ";
  });

  proteinSequence = proteinSequence.replace(/a/g, "ALA");
  proteinSequence = proteinSequence.replace(/b/g, "ASX");
  proteinSequence = proteinSequence.replace(/c/g, "CYS");
  proteinSequence = proteinSequence.replace(/d/g, "ASP");
  proteinSequence = proteinSequence.replace(/e/g, "GLU");
  proteinSequence = proteinSequence.replace(/f/g, "PHE");
  proteinSequence = proteinSequence.replace(/g/g, "GLY");
  proteinSequence = proteinSequence.replace(/h/g, "HIS");
  proteinSequence = proteinSequence.replace(/i/g, "ILE");
  proteinSequence = proteinSequence.replace(/k/g, "LYS");
  proteinSequence = proteinSequence.replace(/l/g, "LEU");
  proteinSequence = proteinSequence.replace(/m/g, "MET");
  proteinSequence = proteinSequence.replace(/n/g, "ASN");
  proteinSequence = proteinSequence.replace(/p/g, "PRO");
  proteinSequence = proteinSequence.replace(/q/g, "GLN");
  proteinSequence = proteinSequence.replace(/r/g, "ARG");
  proteinSequence = proteinSequence.replace(/s/g, "SER");
  proteinSequence = proteinSequence.replace(/t/g, "THR");
  proteinSequence = proteinSequence.replace(/v/g, "VAL");
  proteinSequence = proteinSequence.replace(/w/g, "TRP");
  proteinSequence = proteinSequence.replace(/x/g, "XAA");
  proteinSequence = proteinSequence.replace(/y/g, "TYR");
  proteinSequence = proteinSequence.replace(/z/g, "GLX");
  proteinSequence = proteinSequence.replace(/\*/g, "***");

  proteinSequence = proteinSequence.replace(/\s*(.)(.)(.)\s*/g, function (
    str,
    p1,
    p2,
    p3,
    offset,
    s
  ) {
    return p1 + p2.toLowerCase() + p3.toLowerCase();
  });

  outputWindow.document.write(addReturns(proteinSequence));
  return true;
}
