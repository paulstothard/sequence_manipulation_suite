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

function codonUsage(theDocument) {
  var newDna = "";
  var maxInput = 500000000;
  var codonTable;

  if (testScript() == false) {
    return false;
  }

  if (
    checkFormElement(theDocument.forms[0].elements[0]) == false ||
    checkTextLength(theDocument.forms[0].elements[0].value, maxInput) == false
  ) {
    return false;
  }

  var geneticCode = getGeneticCodeString(
    theDocument.forms[0].elements[4].options[
      theDocument.forms[0].elements[4].selectedIndex
    ].value
  );

  geneticCode = geneticCode.split(/,/);

  if (checkGeneticCode(geneticCode) == false) {
    return false;
  }

  codonTable = makeCodonTable(geneticCode);
  if (codonTable == false) {
    return false;
  }
  resetCounts(codonTable);

  openWindow("Codon Usage");
  openPre();
  var arrayOfFasta = getArrayOfFasta(theDocument.forms[0].elements[0].value);

  for (var i = 0; i < arrayOfFasta.length; i++) {
    newDna = getSequenceFromFasta(arrayOfFasta[i]);
    title = getTitleFromFasta(arrayOfFasta[i]);

    newDna = removeNonDna(newDna);

    outputWindow.document.write(getInfoFromTitleAndSequence(title, newDna));

    addCodons(codonTable, newDna);

    codonTable.determineValues();

    writeCodonTable(codonTable);

    resetCounts(codonTable);

    outputWindow.document.write("\n\n");
  }

  closePre();
  closeWindow();
  return true;
}

function writeCodonTable(codonTable) {
  outputWindow.document.write(
    rightNum("AmAcid", "", 6, "") +
      rightNum("Codon", "", 7, "") +
      rightNum("Number", "", 10, "") +
      rightNum("/1000", "", 12, "") +
      rightNum("Fraction   ..", "", 17, "") +
      "\n\n"
  );
  for (var i = 0; i < codonTable.aminoAcids.length; i++) {
    arrayOfCodons = codonTable.aminoAcids[i].codons;
    for (var j = 0; j < arrayOfCodons.length; j++) {
      outputWindow.document.write(
        rightNum(codonTable.aminoAcids[i].name, "", 3, "") +
          rightNum(arrayOfCodons[j].sequence.toUpperCase(), "", 8, "") +
          rightNum(arrayOfCodons[j].number.toFixed(2), "", 12, "") +
          rightNum(arrayOfCodons[j].perThou.toFixed(2), "", 12, "") +
          rightNum(arrayOfCodons[j].fraction.toFixed(2), "", 12, "") +
          "\n"
      );
    }
    outputWindow.document.write("\n");
  }
}

function addCodons(codonTable, dnaSequence) {
  //replace 'u' with 't'
  dnaSequence = dnaSequence.replace(/u/gi, "t");
  dnaSequence = dnaSequence.replace(/(...)/g, function (str, p1, offset, s) {
    return " " + p1 + " ";
  });

  var matchExp;
  var arrayOfCodons;
  for (var i = 0; i < codonTable.aminoAcids.length; i++) {
    arrayOfCodons = codonTable.aminoAcids[i].codons;
    for (var j = 0; j < arrayOfCodons.length; j++) {
      matchExp = new RegExp(arrayOfCodons[j].sequence, "gi");
      if (dnaSequence.search(matchExp) != -1) {
        arrayOfCodons[j].number =
          arrayOfCodons[j].number + dnaSequence.match(matchExp).length;
      }
    }
  }
}

function resetCounts(codonTable) {
  for (var i = 0; i < codonTable.aminoAcids.length; i++) {
    arrayOfCodons = codonTable.aminoAcids[i].codons;
    for (var j = 0; j < arrayOfCodons.length; j++) {
      arrayOfCodons[j].resetValues();
    }
  }
}

function makeCodonTable(geneticCode) {
  var codonSequence =
    "gggggaggtggcgaggaagatgacgtggtagttgtcgcggcagctgccaggagaagtagcaagaaaaataacatgataattatcacgacaactacctggtgatgttgctagtaatattacttgttatttttctcgtcatcttcccggcgacgtcgccagcaacatcacctgctacttctcccgccacctccc";

  var proteinSequence;
  var codonTable = new CodonTable();

  var geneticCodeMatchExp = getGeneticCodeMatchExp(geneticCode);
  var geneticCodeMatchResult = getGeneticCodeMatchResult(geneticCode);

  codonSequence = codonSequence.replace(/(...)/g, function (
    str,
    p1,
    offset,
    s
  ) {
    return " " + p1 + " ";
  });

  var codonSequenceCopy = codonSequence;

  for (var i = 0; i < geneticCodeMatchExp.length; i++) {
    codonSequence = codonSequence.replace(
      geneticCodeMatchExp[i],
      geneticCodeMatchResult[i]
    );
  }
  var codonArray = codonSequenceCopy.split(/\s+/);

  codonSequence = codonSequence.replace(/\*/g, "Z");
  var proteinArray = codonSequence.split(/\s+/);

  for (var i = 0; i < codonArray.length; i++) {
    //on some systems there will be empty items because of split function
    if (proteinArray[i] == "" && codonArray[i] == "") {
      continue;
    }
    codonTable[proteinArray[i].toLowerCase()].addCodon(
      new Codon(codonArray[i], 0, 0, 0)
    );
  }
  return codonTable;
}

//class CodonTable method determineValues()
function determineValues() {
  var totalAminoAcids = 0;
  var aminoAcid = 0;
  var arrayOfCodons;
  for (var i = 0; i < this.aminoAcids.length; i++) {
    arrayOfCodons = this.aminoAcids[i].codons;
    for (var j = 0; j < arrayOfCodons.length; j++) {
      totalAminoAcids = totalAminoAcids + arrayOfCodons[j].number;
      aminoAcid = aminoAcid + arrayOfCodons[j].number;
    }
    this.aminoAcids[i].count = aminoAcid;
    aminoAcid = 0;
  }

  for (var i = 0; i < this.aminoAcids.length; i++) {
    arrayOfCodons = this.aminoAcids[i].codons;
    for (var j = 0; j < arrayOfCodons.length; j++) {
      if (arrayOfCodons[j].number > 0) {
        arrayOfCodons[j].fraction =
          arrayOfCodons[j].number / this.aminoAcids[i].count;
        arrayOfCodons[j].perThou =
          1000 * (arrayOfCodons[j].number / totalAminoAcids);
      } else {
        arrayOfCodons[j].fraction = 0;
        arrayOfCodons[j].perThou = 0;
      }
    }
  }
}

//class CodonTable
function CodonTable() {
  this.a = new AminoAcid("Ala");
  this.c = new AminoAcid("Cys");
  this.d = new AminoAcid("Asp");
  this.e = new AminoAcid("Glu");
  this.f = new AminoAcid("Phe");
  this.g = new AminoAcid("Gly");
  this.h = new AminoAcid("His");
  this.i = new AminoAcid("Ile");
  this.k = new AminoAcid("Lys");
  this.l = new AminoAcid("Leu");
  this.m = new AminoAcid("Met");
  this.n = new AminoAcid("Asn");
  this.p = new AminoAcid("Pro");
  this.q = new AminoAcid("Gln");
  this.r = new AminoAcid("Arg");
  this.s = new AminoAcid("Ser");
  this.t = new AminoAcid("Thr");
  this.v = new AminoAcid("Val");
  this.w = new AminoAcid("Trp");
  this.y = new AminoAcid("Tyr");
  this.z = new AminoAcid("End");

  this.ala = this.a;
  this.cys = this.c;
  this.asp = this.d;
  this.glu = this.e;
  this.phe = this.f;
  this.gly = this.g;
  this.his = this.h;
  this.ile = this.i;
  this.lys = this.k;
  this.leu = this.l;
  this.met = this.m;
  this.asn = this.n;
  this.pro = this.p;
  this.gln = this.q;
  this.arg = this.r;
  this.ser = this.s;
  this.thr = this.t;
  this.val = this.v;
  this.trp = this.w;
  this.tyr = this.y;
  this.end = this.z;

  this.aminoAcids = new Array(
    this.a,
    this.c,
    this.d,
    this.e,
    this.f,
    this.g,
    this.h,
    this.i,
    this.k,
    this.l,
    this.m,
    this.n,
    this.p,
    this.q,
    this.r,
    this.s,
    this.t,
    this.v,
    this.w,
    this.y,
    this.z
  );
}

//create and throw away a prototype object
new CodonTable();

//define object methods
CodonTable.prototype.determineValues = determineValues;

//class AminoAcid method count()
function count() {
  for (var i = 0; i < this.codons.length; i++) {
    this.number = this.number + this.codons[i];
  }
}

//class AminoAcid method addCodon()
function addCodon(codon) {
  this.codons.push(codon);
}

//class AminoAcid
function AminoAcid(name) {
  this.name = name;
  this.codons = new Array();
  this.number = 0;
}

//create and throw away a prototype object
new AminoAcid("");

//define object methods
AminoAcid.prototype.addCodon = addCodon;
AminoAcid.prototype.count = count;

//class Codon method resetValues()
function resetValues() {
  this.number = 0;
  this.perThou = 0;
  this.fraction = 0;
}

//class Codon
function Codon(sequence, number, perThou, fraction) {
  this.sequence = sequence.toLowerCase();
  this.number = number;
  this.perThou = perThou;
  this.fraction = fraction;
}

//create and throw away a prototype object
new Codon("", 0, 0, 0);

//define object methods
Codon.prototype.resetValues = resetValues;
