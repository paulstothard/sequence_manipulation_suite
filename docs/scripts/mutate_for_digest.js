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

function mutateForDigest(theDocument) {
  var newDna = "";
  var mutatedDna = "";
  var title = "";
  var maxInput = 10000000;
  var TOPOLOGY = "linear";

  if (testScript() == false) {
    return false;
  }

  var restrictionSiteCollection;
  var mutatedRestrictionSitesCollection;

  var geneticCode = getGeneticCodeString(
    theDocument.forms[0].elements[7].options[
      theDocument.forms[0].elements[7].selectedIndex
    ].value
  );

  var restrictionSites =
    theDocument.forms[0].elements[4].options[
      theDocument.forms[0].elements[4].selectedIndex
    ].value;

  if (
    checkFormElement(theDocument.forms[0].elements[0]) == false ||
    checkSequenceLength(theDocument.forms[0].elements[0].value, maxInput) ==
      false
  ) {
    return false;
  }

  geneticCode = geneticCode.split(/,/);

  if (checkGeneticCode(geneticCode) == false) {
    return false;
  }

  restrictionSites = restrictionSites.split(/,/);
  if (checkRestPatterns(restrictionSites) == false) {
    return false;
  }
  var mutatedRestrictionSites = buildMutatedRestrictionSites(restrictionSites);

  openWindow("Mutate for Digest");
  openPre();
  outputWindow.document.write(
    '<span class="mutated_sequence">' +
      "sequence and translations for mutated version" +
      "</span>\n"
  );
  outputWindow.document.write(
    '<span class="current_sequence">' +
      "sequence and translations for non-mutated version" +
      "</span>\n"
  );
  outputWindow.document.write("\n");
  var arrayOfFasta = getArrayOfFasta(theDocument.forms[0].elements[0].value);

  for (var i = 0; i < arrayOfFasta.length; i++) {
    newDna = getSequenceFromFasta(arrayOfFasta[i]);
    title = getTitleFromFasta(arrayOfFasta[i]);

    newDna = removeNonDna(newDna);

    outputWindow.document.write(
      getInfoFromTitleAndSequenceAndTopology(title, newDna, TOPOLOGY)
    );

    restrictionSiteCollection = findRestrictionSites(
      newDna,
      restrictionSites,
      TOPOLOGY
    );
    restrictionSiteCollection.sortRestrictionSites();

    //build restrictionSiteCollection from mutated sites
    mutatedRestrictionSiteCollection = findRestrictionSites(
      newDna,
      mutatedRestrictionSites,
      TOPOLOGY
    );

    mutatedRestrictionSiteCollection = removeNormalMatchesFromMutantSites(
      mutatedRestrictionSiteCollection,
      restrictionSiteCollection
    );

    mutatedRestrictionSiteCollection = removeOverlappingMatchesFromMutantSites(
      mutatedRestrictionSiteCollection
    );

    //use these sites to build mutatedDna
    //modify call below once these are available
    mutatedDna = buildMutatedDna(newDna, mutatedRestrictionSiteCollection);

    layoutRestTrans(
      newDna,
      mutatedDna,
      geneticCode,
      restrictionSiteCollection,
      mutatedRestrictionSiteCollection,
      theDocument.forms[0].elements[5].options[
        theDocument.forms[0].elements[5].selectedIndex
      ].value,
      theDocument.forms[0].elements[6].options[
        theDocument.forms[0].elements[6].selectedIndex
      ].value
    );

    outputWindow.document.write("\n\n");
  }

  closePre();
  closeWindow();
  return true;
}

function layoutRestTrans(
  dnaSequence,
  mutatedDnaSequence,
  geneticCode,
  restrictionSiteCollection,
  mutatedRestrictionSiteCollection,
  basesPerLine,
  readingFrame
) {
  basesPerLine = parseInt(basesPerLine);

  var geneticCodeMatchExp = getGeneticCodeMatchExp(geneticCode);
  var geneticCodeMatchResult = getGeneticCodeMatchResult(geneticCode);
  var spaceString =
    "                                                                                                                                  ";

  var textLayout = new TextLayout();

  if (readingFrame == "3" || readingFrame == "all_three") {
    var translationMut = new TranslationComponent();
    translationMut.spanStart = '<span class="mutated_sequence">';
    translationMut.spanEnd = "</span>";
    translationMut.setCharacters(
      "  " +
        translate(
          mutatedDnaSequence.substring(2, dnaSequence.length),
          geneticCodeMatchExp,
          geneticCodeMatchResult
        )
    );
    textLayout.addComponent(translationMut);

    var translation = new TranslationComponent();
    translation.spanStart = '<span class="current_sequence">';
    translation.spanEnd = "</span>";
    translation.setCharacters(
      "  " +
        translate(
          dnaSequence.substring(2, dnaSequence.length),
          geneticCodeMatchExp,
          geneticCodeMatchResult
        )
    );
    textLayout.addComponent(translation);
  }

  if (readingFrame == "2" || readingFrame == "all_three") {
    var translationMut = new TranslationComponent();
    translationMut.spanStart = '<span class="mutated_sequence">';
    translationMut.spanEnd = "</span>";
    translationMut.setCharacters(
      " " +
        translate(
          mutatedDnaSequence.substring(1, dnaSequence.length),
          geneticCodeMatchExp,
          geneticCodeMatchResult
        )
    );
    textLayout.addComponent(translationMut);

    var translation = new TranslationComponent();
    translation.spanStart = '<span class="current_sequence">';
    translation.spanEnd = "</span>";
    translation.setCharacters(
      " " +
        translate(
          dnaSequence.substring(1, dnaSequence.length),
          geneticCodeMatchExp,
          geneticCodeMatchResult
        )
    );
    textLayout.addComponent(translation);
  }

  if (readingFrame == "1" || readingFrame == "all_three") {
    var translationMut = new TranslationComponent();
    translationMut.spanStart = '<span class="mutated_sequence">';
    translationMut.spanEnd = "</span>";
    translationMut.setCharacters(
      translate(mutatedDnaSequence, geneticCodeMatchExp, geneticCodeMatchResult)
    );
    textLayout.addComponent(translationMut);

    var translation = new TranslationComponent();
    translation.spanStart = '<span class="current_sequence">';
    translation.spanEnd = "</span>";
    translation.setCharacters(
      translate(dnaSequence, geneticCodeMatchExp, geneticCodeMatchResult)
    );
    textLayout.addComponent(translation);
  }

  if (readingFrame == "uppercase") {
    var translationMut = new UppercaseTranslationComponent();
    translationMut.spanStart = '<span class="mutated_sequence">';
    translationMut.spanEnd = "</span>";
    translationMut.setCharacters(
      uppercaseTranslate(
        mutatedDnaSequence,
        geneticCodeMatchExp,
        geneticCodeMatchResult
      )
    );
    textLayout.addComponent(translationMut);

    var translation = new UppercaseTranslationComponent();
    translation.spanStart = '<span class="current_sequence">';
    translation.spanEnd = "</span>";
    translation.setCharacters(
      uppercaseTranslate(
        dnaSequence,
        geneticCodeMatchExp,
        geneticCodeMatchResult
      )
    );
    textLayout.addComponent(translation);
  }

  var dnaMut = new DnaComponent();
  dnaMut.spanStart = '<span class="mutated_sequence">';
  dnaMut.spanEnd = "</span>";
  dnaMut.setCharacters(mutatedDnaSequence);
  textLayout.addComponent(dnaMut);

  var dna = new DnaComponent();
  dna.spanStart = '<span class="current_sequence">';
  dna.spanEnd = "</span>";
  dna.setCharacters(dnaSequence);
  textLayout.addComponent(dna);

  var sitesInRange = new Array();
  var sitesInRangeMut = new Array();

  for (var i = 0; i < dnaSequence.length; i = i + basesPerLine) {
    //first get restriction sites in this range
    sitesInRange = restrictionSiteCollection.getSitesInRange(
      i,
      i + basesPerLine
    );
    sitesInRangeMut = mutatedRestrictionSiteCollection.getSitesInRange(
      i,
      i + basesPerLine
    );
    //reversing gives better appearance
    sitesInRange.reverse();
    sitesInRangeMut.reverse();

    for (var j = 0; j < sitesInRangeMut.length; j++) {
      //9 is added because of indent to dna sequence
      outputWindow.document.write(
        spaceString.substring(0, sitesInRangeMut[j].position - i + 9) +
          '<span class="mutated_sequence">' +
          sitesInRangeMut[j].label +
          "</span>\n"
      );
    }

    for (var j = 0; j < sitesInRange.length; j++) {
      //9 is added because of indent to dna sequence
      outputWindow.document.write(
        spaceString.substring(0, sitesInRange[j].position - i + 9) +
          '<span class="current_sequence">' +
          sitesInRange[j].label +
          "</span>\n"
      );
    }

    textLayout.writeLayout(i, i + basesPerLine);
  }

  return true;
}

function translate(dnaSequence, geneticCodeMatchExp, geneticCodeMatchResult) {
  //don't translate if fewer than three bases
  if (dnaSequence.replace(/[^A-Za-z]/g, "").length < 3) {
    return "";
  }

  dnaSequence = dnaSequence.replace(/(...)/g, function (str, p1, offset, s) {
    return " " + p1 + " ";
  });

  for (var i = 0; i < geneticCodeMatchExp.length; i++) {
    dnaSequence = dnaSequence.replace(
      geneticCodeMatchExp[i],
      geneticCodeMatchResult[i]
    );
  }

  dnaSequence = dnaSequence.replace(/\S{3}/g, "X");
  dnaSequence = dnaSequence.replace(/\s\S{1,2}$/, "");
  dnaSequence = dnaSequence.replace(/^\s/, "");
  return dnaSequence;
}

function uppercaseTranslate(
  dnaSequence,
  geneticCodeMatchExp,
  geneticCodeMatchResult
) {
  dnaSequence = dnaSequence.replace(/[a-z]/g, " ");

  //don't translate if fewer than three bases
  if (dnaSequence.replace(/[^A-Za-z]/g, "").length < 3) {
    return "";
  }

  dnaSequence = dnaSequence.replace(
    /([A-Z])(\s*)([A-Z])(\s*)([A-Z])(\s*)/g,
    function (str, p1, p2, p3, p4, p5, p6, offset, s) {
      return " " + p1 + p3 + p5 + " " + p2 + p4 + p6;
    }
  );

  dnaSequence = dnaSequence.replace(/\s\S{1,2}\s/g, "");

  for (var i = 0; i < geneticCodeMatchExp.length; i++) {
    dnaSequence = dnaSequence.replace(
      geneticCodeMatchExp[i],
      geneticCodeMatchResult[i]
    );
  }

  dnaSequence = dnaSequence.replace(/\S{3}/g, "X");
  dnaSequence = dnaSequence.replace(/^\s/, "");

  return dnaSequence;
}

function buildMutatedRestrictionSites(restrictionSites) {
  var mutatedRestrictionSites = new Array();
  for (var i = 0; i < restrictionSites.length; i++) {
    var site = restrictionSites[i]
      .match(/\/.+\//)
      .toString()
      .replace(/[\/\\]/g, "")
      .toLowerCase();
    var label = restrictionSites[i].match(/\([^\(]+\)/).toString();
    var cutDistance = parseFloat(
      restrictionSites[i]
        .match(/\)\D*\d+/)
        .toString()
        .replace(/\)\D*/, "")
    );
    var singleDegenSites = new Array();
    var doubleDegenSites = new Array();
    for (var j = 0; j < site.length; j++) {
      if (site.charAt(j) != "n" && site.charAt(j) != "N") {
        var pre = site.substring(0, j);
        var post = site.substring(j + 1, site.length);
        var newSite = pre + "N" + post;
        singleDegenSites.push(newSite);
      }
    }

    if (site.length > 6) {
      for (var j = 0; j < singleDegenSites.length; j++) {
        for (var k = 0; k < singleDegenSites[j].length; k++) {
          if (
            singleDegenSites[j].charAt(k) != "n" &&
            singleDegenSites[j].charAt(k) != "N"
          ) {
            var pre = singleDegenSites[j].substring(0, k);
            var post = singleDegenSites[j].substring(
              k + 1,
              singleDegenSites[j].length
            );
            var newSite = pre + "N" + post;
            doubleDegenSites.push(newSite);
          }
        }
      }
    }

    for (var j = 0; j < singleDegenSites.length; j++) {
      mutatedRestrictionSites.push(
        "/" + singleDegenSites[j] + "/" + " " + label + cutDistance
      );
    }

    for (var j = 0; j < doubleDegenSites.length; j++) {
      mutatedRestrictionSites.push(
        "/" + doubleDegenSites[j] + "/" + " " + label + cutDistance
      );
    }
  }
  return mutatedRestrictionSites;
}

function removeNormalMatchesFromMutantSites(
  mutatedRestrictionSiteCollection,
  restrictionSiteCollection
) {
  var originalMutatedRestrictionSites = new Array();
  for (
    var i = 0;
    i < mutatedRestrictionSiteCollection.restrictionSites.length;
    i++
  ) {
    var mutatedSite = mutatedRestrictionSiteCollection.restrictionSites[i];
    var isOriginal = true;
    for (
      var j = 0;
      j < restrictionSiteCollection.restrictionSites.length;
      j++
    ) {
      var normalSite = restrictionSiteCollection.restrictionSites[j];
      if (normalSite.position == mutatedSite.position) {
        isOriginal = false;
        break;
      }
    }
    if (isOriginal) {
      originalMutatedRestrictionSites.push(mutatedSite);
    }
  }
  mutatedRestrictionSiteCollection.restrictionSites = originalMutatedRestrictionSites;
  return mutatedRestrictionSiteCollection;
}

function removeOverlappingMatchesFromMutantSites(
  mutatedRestrictionSiteCollection
) {
  var originalMutatedRestrictionSites = new Array();
  for (
    var i = 0;
    i < mutatedRestrictionSiteCollection.restrictionSites.length;
    i++
  ) {
    var mutatedSite = mutatedRestrictionSiteCollection.restrictionSites[i];
    var isOriginal = true;
    for (var j = 0; j < originalMutatedRestrictionSites.length; j++) {
      var mutatedSiteJ = originalMutatedRestrictionSites[j];
      //don't want one site to alter protein translation under another site
      var startRangeJ =
        mutatedSiteJ.position +
        mutatedSiteJ.cutDistance -
        mutatedSiteJ.iupacPattern.length -
        2;
      var endRangeJ = mutatedSiteJ.position + mutatedSiteJ.cutDistance - 1 + 2;

      var startRange =
        mutatedSite.position +
        mutatedSite.cutDistance -
        mutatedSite.iupacPattern.length;
      var endRange = mutatedSite.position + mutatedSite.cutDistance - 1;

      if (startRange <= startRangeJ && endRange >= startRangeJ) {
        isOriginal = false;
        break;
      }

      if (startRange <= endRangeJ && endRange >= endRangeJ) {
        isOriginal = false;
        break;
      }

      if (startRange <= startRangeJ && endRange >= endRangeJ) {
        isOriginal = false;
        break;
      }

      if (startRange >= startRangeJ && endRange <= endRangeJ) {
        isOriginal = false;
        break;
      }
    }
    if (isOriginal) {
      originalMutatedRestrictionSites.push(mutatedSite);
    }
  }
  mutatedRestrictionSiteCollection.restrictionSites = originalMutatedRestrictionSites;
  return mutatedRestrictionSiteCollection;
}

function buildMutatedDna(originalDna, mutatedRestrictionSiteCollection) {
  var mutatedDna = originalDna;
  var mutatedDnaArray = new Array();
  mutatedRestrictionSiteCollection.sortRestrictionSites();
  mutatedRestrictionSiteCollection.restrictionSites.reverse();
  for (
    var i = 0;
    i < mutatedRestrictionSiteCollection.restrictionSites.length;
    i++
  ) {
    var mutatedSite = mutatedRestrictionSiteCollection.restrictionSites[i];
    var siteStart =
      mutatedSite.position +
      mutatedSite.cutDistance -
      mutatedSite.iupacPattern.length;
    var siteEnd = mutatedSite.position + mutatedSite.cutDistance - 1;
    var siteLength = siteEnd - siteStart;

    mutatedDnaArray.push(
      mutatedDna.substring(0, siteStart) +
        replaceMutatedDnaSegment(
          mutatedDna.substring(siteStart, siteEnd + 1),
          mutatedSite.iupacPattern,
          mutatedSite.label
        )
    );
    mutatedDnaArray.push(mutatedDna.substring(siteEnd + 1, mutatedDna.length));

    mutatedDna = mutatedDnaArray.join("");
    mutatedDnaArray = new Array();
  }
  mutatedRestrictionSiteCollection.restrictionSites.reverse();
  return mutatedDna;
}

function replaceMutatedDnaSegment(originalSegment, iupacPattern, label) {
  var newSegment = "";
  var genericSite;
  //determine genericSite from label
  genericSite = label
    .match(/[a-z\|]+\s\d+/)
    .toString()
    .replace(/\||\s\d+/g, "");
  var isUpperCase;
  var charToAdd;
  for (var i = 0; i < originalSegment.length; i++) {
    if (originalSegment.charAt(i).toString().search(/[A-Z]/) == -1) {
      isUpperCase = false;
    } else {
      isUpperCase = true;
    }

    if (originalSegment.charAt(i) == genericSite.charAt(i)) {
      charToAdd = originalSegment.charAt(i);
    } else if (iupacPattern.charAt(i) == "N") {
      charToAdd = genericSite.charAt(i);
    } else {
      charToAdd = originalSegment.charAt(i);
    }

    if (isUpperCase) {
      newSegment = newSegment + charToAdd.toString().toUpperCase();
    } else {
      newSegment = newSegment + charToAdd.toString().toLowerCase();
    }
  }

  return newSegment;
}

function findRestrictionSites(sequence, arrayOfItems, dnaConformation) {
  var lookAhead = 50;
  var lowerLimit = 0;
  var upperLimit = sequence.length;
  var shiftValue = 0;
  var cutDistance;
  var matchExp;
  var matchPosition;
  var matchArray;
  var label;
  var timesFound = 0;
  var tempArray = new Array();

  var restrictionSiteCollection = new RestrictionSiteCollection();

  if (dnaConformation == "circular") {
    shiftValue = sequence.substring(0, lookAhead).length;
    sequence =
      sequence.substring(sequence.length - lookAhead, sequence.length) +
      sequence +
      sequence.substring(0, lookAhead);
    lowerLimit = 0 + shiftValue;
    upperLimit = upperLimit + shiftValue;
  }

  for (var i = 0; i < arrayOfItems.length; i++) {
    var iupacPattern = arrayOfItems[i]
      .match(/\/.+\//)
      .toString()
      .replace(/[\/\\]/g, "");
    matchExp = "/" + convertDegenerates(iupacPattern) + "/" + "gi";
    matchPosition = 0;
    matchExp = eval(matchExp);
    cutDistance = parseFloat(
      arrayOfItems[i]
        .match(/\)\D*\d+/)
        .toString()
        .replace(/\)\D*/, "")
    );
    label = arrayOfItems[i]
      .match(/\([^\(]+\)/)
      .toString()
      .replace(/\(|\)/g, "");

    while ((matchArray = matchExp.exec(sequence))) {
      matchPosition = matchExp.lastIndex;
      matchPosition = matchPosition - cutDistance;
      if (matchPosition >= lowerLimit && matchPosition < upperLimit) {
        timesFound++;
        tempArray.push(
          new RestrictionSite(
            label + " " + (matchPosition - shiftValue + 1),
            matchPosition - shiftValue,
            cutDistance,
            iupacPattern
          )
        );
      }

      matchExp.lastIndex = matchExp.lastIndex - RegExp.lastMatch.length + 1;
    }

    for (var j = 0; j < tempArray.length; j++) {
      tempArray[j].numberOfCuts = timesFound;
      restrictionSiteCollection.addRestrictionSite(tempArray[j]);
    }
    timesFound = 0;
    tempArray = new Array();
  }

  return restrictionSiteCollection;
}

//------------------------------------ TextLayout Class
//TextLayout writeLayout method
function writeLayout(start, stop) {
  for (var i = 0; i < this.components.length; i++) {
    this.components[i].writeLayoutComponent(start, stop);
  }
}

//TextLayout addComponent method
function addComponent(component) {
  this.components.push(component);
}

//TextLayout class
function TextLayout() {
  this.components = new Array();
}

//create and throw away a prototype object
new TextLayout();

//define object methods
TextLayout.prototype.writeLayout = writeLayout;
TextLayout.prototype.addComponent = addComponent;
//------------------------------------

//------------------------------------ LayoutComponent Abstract Class
//LayoutComponent writeLayoutComponent method
function writeLayoutComponent(start, stop) {
  //abstract method
}

//LayoutComponent setCharacters method
function setCharacters(text) {
  if (text.search(/./) != -1) {
    this.characters = text.match(/./g);
  }
}

//LayoutComponent isRoom method
function isRoom(start, stop) {
  //return true if nothing in characters array in this range
  var rangeToCheck = this.characters.slice(start, stop);
  rangeToCheck = rangeToCheck.join("");
  if (rangeToCheck.search(/\w/) == -1) {
    return true;
  } else {
    return false;
  }
}

//LayoutComponent class
function LayoutComponent() {
  this.characters = new Array();
  this.positionLabel = 1;
  this.spanStart = "";
  this.spanEnd = "";
}

//create and throw away a prototype object
new LayoutComponent();

//define object methods
LayoutComponent.prototype.writeLayoutComponent = writeLayoutComponent;
LayoutComponent.prototype.setCharacters = setCharacters;
LayoutComponent.prototype.isRoom = isRoom;
//------------------------------------

//------------------------------------ TranslationComponent Class extend LayoutComponent Class
//constructor
function TranslationComponent() {
  this.characters = new Array();
  this.positionLabel = 1;
  this.spanStart = "";
  this.spanEnd = "";
}

TranslationComponent.prototype = new LayoutComponent();

//override writeLayoutComponent
TranslationComponent.prototype.writeLayoutComponent = function (start, stop) {
  outputWindow.document.write(rightNum(this.positionLabel, "", 8, ""));
  outputWindow.document.write(this.spanStart);
  outputWindow.document.write(this.characters.slice(start, stop).join(""));
  outputWindow.document.write(this.spanEnd + "\n");
  this.positionLabel = this.positionLabel + (stop - start) / 3;
};

//------------------------------------

//------------------------------------ UppercaseTranslationComponent Class extend LayoutComponent Class
//constructor
function UppercaseTranslationComponent() {
  this.characters = new Array();
  this.positionLabel = 1;
  this.spanStart = "";
  this.spanEnd = "";
}

UppercaseTranslationComponent.prototype = new LayoutComponent();

//override writeLayoutComponent
UppercaseTranslationComponent.prototype.writeLayoutComponent = function (
  start,
  stop
) {
  var textToWrite = this.characters.slice(start, stop).join("");
  if (textToWrite.search(/\w/) != -1) {
    outputWindow.document.write(rightNum(this.positionLabel, "", 8, ""));
    outputWindow.document.write(this.spanStart);
    this.positionLabel =
      this.positionLabel + textToWrite.match(/[A-Z]/g).length;
    outputWindow.document.write(textToWrite);
    outputWindow.document.write(this.spanEnd + "\n");
  }
};
//------------------------------------

//------------------------------------ DnaComponent Class extend LayoutComponent Class
//constructor
function DnaComponent() {
  this.characters = new Array();
  this.positionLabel = 1;
  this.spanStart = "";
  this.spanEnd = "";
}
DnaComponent.prototype = new LayoutComponent();

//override writeLayoutComponent
DnaComponent.prototype.writeLayoutComponent = function (start, stop) {
  outputWindow.document.write(rightNum(this.positionLabel, "", 8, ""));
  outputWindow.document.write(this.spanStart);
  outputWindow.document.write(this.characters.slice(start, stop).join(""));
  outputWindow.document.write(this.spanEnd + "\n");
  this.positionLabel = this.positionLabel + (stop - start);
};
//------------------------------------

//------------------------------------ RestrictionSite class
//RestrictionSite class
function RestrictionSite(label, position, cutDistance, iupacPattern) {
  this.label = label;
  this.position = position;
  this.cutDistance = cutDistance;
  this.iupacPattern = iupacPattern;
  this.numberOfCuts;
}
//------------------------------------

//------------------------------------ RestrictionSiteCollection class
//RestrictionSiteCollection addRestrictionSite method
function addRestrictionSite(restrictionSite) {
  this.restrictionSites.push(restrictionSite);
}

//RestrictionSiteCollection sortRestrictionSites method
function sortRestrictionSites() {
  this.restrictionSites.sort(restrictionSiteSorter);
}

//RestrictionSiteCollection getSitesInRange method
function getSitesInRange(start, stop) {
  var arrayToReturn = new Array();
  //start at end
  for (var i = this.restrictionSites.length - 1; i >= 0; i--) {
    if (
      this.restrictionSites[i].position >= start &&
      this.restrictionSites[i].position < stop
    ) {
      arrayToReturn.push(this.restrictionSites.pop());
    } else {
      break;
    }
  }
  return arrayToReturn;
}

//RestrictionSiteCollection class
function RestrictionSiteCollection() {
  this.restrictionSites = new Array();
}

//create and throw away a prototype object
new RestrictionSiteCollection();

//define object methods
RestrictionSiteCollection.prototype.addRestrictionSite = addRestrictionSite;
RestrictionSiteCollection.prototype.sortRestrictionSites = sortRestrictionSites;
RestrictionSiteCollection.prototype.getSitesInRange = getSitesInRange;

//------------------------------------

//sort so first ones at end
function restrictionSiteSorter(a, b) {
  if (a.position < b.position) {
    return 1;
  }
  if (a.position > b.position) {
    return -1;
  } else {
    return 0;
  }
}
