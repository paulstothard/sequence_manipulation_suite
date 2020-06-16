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

function transMap(theDocument) {
  var newDna = "";
  var title = "";
  var maxInput = 500000000;

  if (testScript() == false) {
    return false;
  }

  var geneticCode = getGeneticCodeString(
    theDocument.forms[0].elements[6].options[
      theDocument.forms[0].elements[6].selectedIndex
    ].value
  );

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

  openWindow("Translation Map");
  openPre();
  var arrayOfFasta = getArrayOfFasta(theDocument.forms[0].elements[0].value);

  for (var i = 0; i < arrayOfFasta.length; i++) {
    newDna = getSequenceFromFasta(arrayOfFasta[i]);
    title = getTitleFromFasta(arrayOfFasta[i]);

    newDna = removeNonDna(newDna);

    outputWindow.document.write(getInfoFromTitleAndSequence(title, newDna));

    layoutTranslation(
      newDna,
      geneticCode,
      theDocument.forms[0].elements[4].options[
        theDocument.forms[0].elements[4].selectedIndex
      ].value,
      theDocument.forms[0].elements[5].options[
        theDocument.forms[0].elements[5].selectedIndex
      ].value
    );

    outputWindow.document.write("\n\n");
  }

  closePre();
  closeWindow();
  return true;
}

function layoutTranslation(
  dnaSequence,
  geneticCode,
  basesPerLine,
  readingFrame
) {
  basesPerLine = parseInt(basesPerLine);

  var geneticCodeMatchExp = getGeneticCodeMatchExp(geneticCode);
  var geneticCodeMatchResult = getGeneticCodeMatchResult(geneticCode);

  var textLayout = new TextLayout();

  if (readingFrame == "3" || readingFrame == "all_three") {
    var translation = new TranslationComponent();
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
    var translation = new TranslationComponent();
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
    var translation = new TranslationComponent();
    translation.setCharacters(
      translate(dnaSequence, geneticCodeMatchExp, geneticCodeMatchResult)
    );
    textLayout.addComponent(translation);
  }

  if (readingFrame == "uppercase") {
    var translation = new UppercaseTranslationComponent();
    translation.setCharacters(
      uppercaseTranslate(
        dnaSequence,
        geneticCodeMatchExp,
        geneticCodeMatchResult
      )
    );
    textLayout.addComponent(translation);
  }

  var dna = new DnaComponent();
  dna.setCharacters(dnaSequence);
  textLayout.addComponent(dna);

  var ruler = new RulerComponent();
  ruler.setCharacters(dnaSequence);
  ruler.buildRuler();
  textLayout.addComponent(ruler);

  dna = new DnaComponent();
  dna.setCharacters(complement(dnaSequence));
  textLayout.addComponent(dna);

  for (var i = 0; i < dnaSequence.length; i = i + basesPerLine) {
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
}

TranslationComponent.prototype = new LayoutComponent();

//override writeLayoutComponent
TranslationComponent.prototype.writeLayoutComponent = function (start, stop) {
  outputWindow.document.write(rightNum(this.positionLabel, "", 8, ""));
  outputWindow.document.write(
    this.characters.slice(start, stop).join("") + "\n"
  );
  this.positionLabel = this.positionLabel + (stop - start) / 3;
};

//------------------------------------

//------------------------------------ UppercaseTranslationComponent Class extend LayoutComponent Class
//constructor
function UppercaseTranslationComponent() {
  this.characters = new Array();
  this.positionLabel = 1;
}

UppercaseTranslationComponent.prototype = new LayoutComponent();

//override writeLayoutComponent
UppercaseTranslationComponent.prototype.writeLayoutComponent = function (
  start,
  stop
) {
  var textToWrite = this.characters.slice(start, stop).join("") + "\n";
  if (textToWrite.search(/\w/) != -1) {
    outputWindow.document.write(rightNum(this.positionLabel, "", 8, ""));
    this.positionLabel =
      this.positionLabel + textToWrite.match(/[A-Z]/g).length;
    outputWindow.document.write(textToWrite);
  }
};
//------------------------------------

//------------------------------------ DnaComponent Class extend LayoutComponent Class
//constructor
function DnaComponent() {
  this.characters = new Array();
  this.positionLabel = 1;
}
DnaComponent.prototype = new LayoutComponent();

//override writeLayoutComponent
DnaComponent.prototype.writeLayoutComponent = function (start, stop) {
  outputWindow.document.write(rightNum(this.positionLabel, "", 8, ""));
  outputWindow.document.write(
    this.characters.slice(start, stop).join("") + "\n"
  );
  this.positionLabel = this.positionLabel + (stop - start);
};
//------------------------------------

//------------------------------------ RulerComponent Class extend LayoutComponent Class
//constructor
function RulerComponent() {
  this.characters = new Array();
  this.positionLabel = 1;
}
RulerComponent.prototype = new LayoutComponent();

//override writeLayoutComponent
RulerComponent.prototype.writeLayoutComponent = function (start, stop) {
  outputWindow.document.write(rightNum(this.positionLabel, "", 8, ""));
  var text = this.characters.slice(start, stop).join("");
  text = text.replace(/^(\d+)/g, function (str, p1, offset, s) {
    return p1.replace(/./g, " ");
  });
  text = text.replace(/(\d+)$/g, function (str, p1, offset, s) {
    return p1.replace(/./g, " ");
  });
  outputWindow.document.write(text + "\n");
  this.positionLabel = this.positionLabel + (stop - start);
};

//RulerComponent buildRuler method
function buildRuler() {
  //do something
  var sequence = this.characters.join("");
  var count = 0;
  var spacing = "         ";
  sequence = sequence.replace(/(.{1,10})/g, function (str, p1, offset, s) {
    var ruler = count + spacing;
    if (count == 0) {
      ruler = spacing;
    }
    count = count + 10;
    return ruler.substring(0, 10);
  });
  this.characters = sequence.match(/./g);
}

//create and throw away a prototype object
new RulerComponent();

//define object methods
RulerComponent.prototype.buildRuler = buildRuler;

//------------------------------------
