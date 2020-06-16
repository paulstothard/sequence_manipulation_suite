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

function restSummary(theDocument) {
  var newDna = "";
  var title = "";
  var maxInput = 100000000;

  var restrictionSites = getRestrictionSiteString("standard");

  if (
    checkFormElement(theDocument.forms[0].elements[0]) == false ||
    checkSequenceLength(theDocument.forms[0].elements[0].value, maxInput) ==
      false
  ) {
    return false;
  }

  itemsToCheck = restrictionSites.split(/,/);
  if (checkRestPatterns(itemsToCheck) == false) {
    return false;
  }

  openWindow("Restriction Summary");
  //openPre();
  outputWindow.document.write(
    '<span class="one">' + "cuts once" + "</span><br />\n"
  );
  outputWindow.document.write(
    '<span class="two">' + "cuts twice" + "</span><br />\n"
  );
  outputWindow.document.write("\n");
  var arrayOfFasta = getArrayOfFasta(theDocument.forms[0].elements[0].value);

  for (var i = 0; i < arrayOfFasta.length; i++) {
    newDna = getSequenceFromFasta(arrayOfFasta[i]);
    title = getTitleFromFasta(arrayOfFasta[i]);

    newDna = removeNonDna(newDna);

    outputWindow.document.write(
      getInfoFromTitleAndSequenceAndTopology(
        title,
        newDna,
        theDocument.forms[0].elements[4].options[
          theDocument.forms[0].elements[4].selectedIndex
        ].value
      )
    );

    writeRestrictionSites(
      newDna,
      itemsToCheck,
      theDocument.forms[0].elements[4].options[
        theDocument.forms[0].elements[4].selectedIndex
      ].value
    );

    outputWindow.document.write("<br />\n<br />\n");
  }

  //closePre();
  closeWindow();
  return true;
}
