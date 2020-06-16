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

function identSim(theDocument) {
  var maxInput = 20000000;
  var theAlignment = "";
  var alignArray = new Array();
  var groupString = "";
  var arrayOfGroups = new Array();

  var titleArray = new Array();
  var sequenceArray = new Array();

  var longestTitle;

  if (testScript() == false) {
    return false;
  }

  if (
    checkFormElement(theDocument.forms[0].elements[0]) == false ||
    checkTextLength(theDocument.forms[0].elements[0].value, maxInput) == false
  ) {
    return false;
  }

  theAlignment = "X" + theDocument.forms[0].elements[0].value;
  alignArray = theAlignment.split(/[>%#]/);

  if (earlyCheckAlign(alignArray) == false) {
    return false;
  }

  for (var i = 1; i < alignArray.length; i++) {
    titleArray[i - 1] = alignArray[i].match(/[^\f\n\r]+[\f\n\r]/);
    titleArray[i - 1] = filterFastaTitle(titleArray[i - 1].toString()).replace(
      /[\f\n\r]/g,
      ""
    );
    titleArray[i - 1] = titleArray[i - 1].substring(0, 20);
    if (i == 1) {
      longestTitle = titleArray[i - 1].length;
    } else if (titleArray[i - 1].length > longestTitle) {
      longestTitle = titleArray[i - 1].length;
    }
    sequenceArray[i - 1] = alignArray[i].replace(/[^\f\n\r]+[\f\n\r]/, "");
    sequenceArray[i - 1] = filterAlignSeq(sequenceArray[i - 1]);
  }

  if (checkAlign(titleArray, sequenceArray) == false) {
    return false;
  }

  groupString = theDocument.forms[0].elements[1].value
    .replace(/\s/g, "")
    .toUpperCase();
  arrayOfGroups = groupString.split(/,/);
  if (checkGroupInput(arrayOfGroups) == false) {
    return false;
  }

  openWindowAlign("Ident and Sim");
  openPre();
  writeIdentAndSim(titleArray, sequenceArray, arrayOfGroups);
  closePre();
  closeWindow();
  return true;
}

function writeIdentAndSim(titleArray, sequenceArray, arrayOfGroups) {
  var identical = 0;
  var similar = 0;
  var alignLength = 0;
  for (var k = 0; k < sequenceArray.length; k++) {
    for (var m = k + 1; m < sequenceArray.length; m++) {
      for (var i = 0; i < sequenceArray[0].length; i++) {
        alignLength = alignLength + 1;
        if (
          sequenceArray[k].charAt(i).toUpperCase() ==
            sequenceArray[m].charAt(i).toUpperCase() &&
          sequenceArray[k].charAt(i).toUpperCase() != "X"
        ) {
          if (
            sequenceArray[k].charAt(i) != "-" &&
            sequenceArray[k].charAt(i) != "."
          ) {
            identical = identical + 1;
          } else {
            alignLength = alignLength - 1;
          }
        } else {
          for (var j = 0; j < arrayOfGroups.length; j++) {
            if (
              arrayOfGroups[j].search(
                sequenceArray[k].charAt(i).toUpperCase()
              ) != -1 &&
              arrayOfGroups[j].search(
                sequenceArray[m].charAt(i).toUpperCase()
              ) != -1
            ) {
              similar = similar + 1;
              break;
            }
          }
        }
      }
      outputWindow.document.write(
        "<b>Results for " + titleArray[k] + " vs " + titleArray[m] + ":</b>\n"
      );
      outputWindow.document.write("  Alignment length: " + alignLength + "\n");
      outputWindow.document.write("Identical residues: " + identical + "\n");
      outputWindow.document.write("  Similar residues: " + similar + "\n");
      if (identical == 0) {
        outputWindow.document.write("  Percent identity: " + 0 + "\n");
      } else {
        outputWindow.document.write(
          "  Percent identity: " +
            ((identical / alignLength) * 100).toFixed(2) +
            "\n"
        );
      }
      if (similar == 0 && identical == 0) {
        outputWindow.document.write("Percent similarity: " + 0 + "\n");
      } else {
        outputWindow.document.write(
          "Percent similarity: " +
            (((identical + similar) / alignLength) * 100).toFixed(2) +
            "\n"
        );
      }
      outputWindow.document.write("\n");
      identical = 0;
      similar = 0;
      alignLength = 0;
    }
  }
  return true;
}
