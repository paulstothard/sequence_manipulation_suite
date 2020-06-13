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

function revTrans(theDocument) {
  var newProtein = "";
  var maxInput = 20000000;

  if (testScript() == false) {
    return false;
  }

  var codonTable;

  if (
    checkFormElement(theDocument.forms[0].elements[0]) == false ||
    checkSequenceLength(theDocument.forms[0].elements[0].value, maxInput) ==
      false ||
    checkCodonTable(theDocument.forms[0].elements[4].value) == false
  ) {
    return false;
  }

  codonTable = makeCodonTable(theDocument.forms[0].elements[4].value);
  if (codonTable == false) {
    return false;
  }

  openWindow("Reverse Translate");
  openPre();
  var arrayOfFasta = getArrayOfFasta(theDocument.forms[0].elements[0].value);

  for (var i = 0; i < arrayOfFasta.length; i++) {
    newProtein = getSequenceFromFasta(arrayOfFasta[i]);
    title = getTitleFromFasta(arrayOfFasta[i]);

    newProtein = removeNonProteinAllowX(newProtein);

    outputWindow.document.write(getInfoFromTitleAndSequence(title, newProtein));
    writeRevTransSeqNoDegen(newProtein, title, codonTable);
    outputWindow.document.write("\n");
    writeRevTransSeqDegen(newProtein, title, codonTable);
    outputWindow.document.write("\n");
    outputWindow.document.write("Graph of base probabilities:\n");
    writeRevTransGraph(newProtein, codonTable);

    outputWindow.document.write("\n\n");
  }

  closePre();
  closeWindow();
  return true;
}

function writeRevTransSeqNoDegen(protein, title, codonTable) {
  var aminoAcid;

  protein = protein.replace(/\*/g, "z");

  protein = protein.replace(/(.)/g, function (str, p1, offset, s) {
    aminoAcid = codonTable[p1.toString().toLowerCase()];
    return aminoAcid.mostCommonCodon;
  });
  outputWindow.document.write(
    "&gt;" +
      "reverse translation of " +
      title +
      " to a " +
      protein.length +
      " base sequence of most likely codons.\n"
  );
  outputWindow.document.write(addReturns(protein));
  outputWindow.document.write("\n");
  return true;
}

function writeRevTransSeqDegen(protein, title, codonTable) {
  var aminoAcid;

  protein = protein.replace(/\*/g, "z");

  protein = protein.replace(/(.)/g, function (str, p1, offset, s) {
    aminoAcid = codonTable[p1.toString().toLowerCase()];
    return aminoAcid.degenCodon;
  });
  outputWindow.document.write(
    "&gt;" +
      "reverse translation of " +
      title +
      " to a " +
      protein.length +
      " base sequence of consensus codons.\n"
  );
  outputWindow.document.write(addReturns(protein));
  outputWindow.document.write("\n");
  return true;
}

function writeRevTransGraph(protein, codonTable) {
  var aminoAcid;

  protein = protein.replace(/\*/g, "z");

  protein = protein.replace(/(.)/g, function (str, p1, offset, s) {
    aminoAcid = codonTable[p1.toString().toLowerCase()];
    return (
      "<b>" +
      (offset + 1) +
      "_" +
      str +
      "_" +
      "first</b>\n" +
      aminoAcid.rulerPosOne +
      "<b>" +
      (offset + 1) +
      "_" +
      str +
      "_" +
      "second</b>\n" +
      aminoAcid.rulerPosTwo +
      "<b>" +
      (offset + 1) +
      "_" +
      str +
      "_" +
      "third</b>\n" +
      aminoAcid.rulerPosThree +
      "\n"
    );
  });

  outputWindow.document.write(protein);

  return true;
}

function makeCodonTable(gcgTable) {
  gcgTable = gcgTable.replace(/[^\.]*\.\./, "");
  var tableArray = gcgTable.split(/[\f\n\r]/);
  var re = /(\w+)\s+(\w+)\s+(\S+)\s+(\S+)\s+(\S+)/g;
  var matchArray;
  var codonTable = new CodonTable();

  for (var i = 0; i < tableArray.length; i++) {
    while ((matchArray = re.exec(tableArray[i]))) {
      try {
        codonTable[matchArray[1].toLowerCase()].addCodon(
          new Codon(
            matchArray[2],
            parseFloat(matchArray[3]),
            parseFloat(matchArray[4]),
            parseFloat(matchArray[5])
          )
        );
      } catch (e) {
        alert(
          "There is a problem with a line of the codon table: " +
            matchArray[1] +
            " " +
            matchArray[2] +
            " " +
            matchArray[3] +
            " " +
            matchArray[4] +
            " " +
            matchArray[4]
        );
        return false;
      }
    }
  }

  codonTable.a.determineBaseFreq();
  codonTable.c.determineBaseFreq();
  codonTable.d.determineBaseFreq();
  codonTable.e.determineBaseFreq();
  codonTable.f.determineBaseFreq();
  codonTable.g.determineBaseFreq();
  codonTable.h.determineBaseFreq();
  codonTable.i.determineBaseFreq();
  codonTable.k.determineBaseFreq();
  codonTable.l.determineBaseFreq();
  codonTable.m.determineBaseFreq();
  codonTable.n.determineBaseFreq();
  codonTable.p.determineBaseFreq();
  codonTable.q.determineBaseFreq();
  codonTable.r.determineBaseFreq();
  codonTable.s.determineBaseFreq();
  codonTable.t.determineBaseFreq();
  codonTable.v.determineBaseFreq();
  codonTable.w.determineBaseFreq();
  codonTable.y.determineBaseFreq();
  codonTable.z.determineBaseFreq();

  //assign equal frequency to bases coding for unknown amino acid
  codonTable.x.determineBaseFreq();
  codonTable.x.baseFreqPosOne = new Array(0.25, 0.25, 0.25, 0.25);
  codonTable.x.baseFreqPosTwo = new Array(0.25, 0.25, 0.25, 0.25);
  codonTable.x.baseFreqPosThree = new Array(0.25, 0.25, 0.25, 0.25);

  codonTable.a.fillRuler();
  codonTable.c.fillRuler();
  codonTable.d.fillRuler();
  codonTable.e.fillRuler();
  codonTable.f.fillRuler();
  codonTable.g.fillRuler();
  codonTable.h.fillRuler();
  codonTable.i.fillRuler();
  codonTable.k.fillRuler();
  codonTable.l.fillRuler();
  codonTable.m.fillRuler();
  codonTable.n.fillRuler();
  codonTable.p.fillRuler();
  codonTable.q.fillRuler();
  codonTable.r.fillRuler();
  codonTable.s.fillRuler();
  codonTable.t.fillRuler();
  codonTable.v.fillRuler();
  codonTable.w.fillRuler();
  codonTable.y.fillRuler();
  codonTable.z.fillRuler();
  codonTable.x.fillRuler();

  return codonTable;
}

//class CodonTable
function CodonTable() {
  this.a = new AminoAcid();
  this.c = new AminoAcid();
  this.d = new AminoAcid();
  this.e = new AminoAcid();
  this.f = new AminoAcid();
  this.g = new AminoAcid();
  this.h = new AminoAcid();
  this.i = new AminoAcid();
  this.k = new AminoAcid();
  this.l = new AminoAcid();
  this.m = new AminoAcid();
  this.n = new AminoAcid();
  this.p = new AminoAcid();
  this.q = new AminoAcid();
  this.r = new AminoAcid();
  this.s = new AminoAcid();
  this.t = new AminoAcid();
  this.v = new AminoAcid();
  this.w = new AminoAcid();
  this.y = new AminoAcid();
  this.z = new AminoAcid();
  this.x = new AminoAcid();

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
  this.xaa = this.x;
}

//class AminoAcid method fillRuler()
function fillRuler() {
  var markG =
    "gggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg";
  var markA =
    "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
  var markT =
    "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
  var markC =
    "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";

  var markLength = markG.length;

  this.rulerPosOne =
    this.getRuler(this.baseFreqPosOne[0], markG) +
    this.getRuler(this.baseFreqPosOne[1], markA) +
    this.getRuler(this.baseFreqPosOne[2], markT) +
    this.getRuler(this.baseFreqPosOne[3], markC);

  this.rulerPosTwo =
    this.getRuler(this.baseFreqPosTwo[0], markG) +
    this.getRuler(this.baseFreqPosTwo[1], markA) +
    this.getRuler(this.baseFreqPosTwo[2], markT) +
    this.getRuler(this.baseFreqPosTwo[3], markC);

  this.rulerPosThree =
    this.getRuler(this.baseFreqPosThree[0], markG) +
    this.getRuler(this.baseFreqPosThree[1], markA) +
    this.getRuler(this.baseFreqPosThree[2], markT) +
    this.getRuler(this.baseFreqPosThree[3], markC);

  //calculate some other values
  this.setMostCommonCodon();
  this.setDegenCodon();
}

//class AminoAcid method getRuler()
function getRuler(freq, mark) {
  return (
    mark.substring(0, 1) +
    mark.substring(0, freq * mark.length) +
    " " +
    freq.toFixed(2) +
    "\n"
  );
}

//class AminoAcid method addCodon()
function addCodon(codon) {
  this.codons.push(codon);
}

//class AminoAcid method determineBaseFreq()
function determineBaseFreq() {
  this.fixFraction();
  for (var i = 0; i < this.codons.length; i++) {
    if (this.codons[i].sequence.charAt(0) == "g") {
      this.baseFreqPosOne[0] = this.baseFreqPosOne[0] + this.codons[i].fraction;
    } else if (this.codons[i].sequence.charAt(0) == "a") {
      this.baseFreqPosOne[1] = this.baseFreqPosOne[1] + this.codons[i].fraction;
    } else if (this.codons[i].sequence.charAt(0) == "t") {
      this.baseFreqPosOne[2] = this.baseFreqPosOne[2] + this.codons[i].fraction;
    } else if (this.codons[i].sequence.charAt(0) == "c") {
      this.baseFreqPosOne[3] = this.baseFreqPosOne[3] + this.codons[i].fraction;
    }

    if (this.codons[i].sequence.charAt(1) == "g") {
      this.baseFreqPosTwo[0] = this.baseFreqPosTwo[0] + this.codons[i].fraction;
    } else if (this.codons[i].sequence.charAt(1) == "a") {
      this.baseFreqPosTwo[1] = this.baseFreqPosTwo[1] + this.codons[i].fraction;
    } else if (this.codons[i].sequence.charAt(1) == "t") {
      this.baseFreqPosTwo[2] = this.baseFreqPosTwo[2] + this.codons[i].fraction;
    } else if (this.codons[i].sequence.charAt(1) == "c") {
      this.baseFreqPosTwo[3] = this.baseFreqPosTwo[3] + this.codons[i].fraction;
    }

    if (this.codons[i].sequence.charAt(2) == "g") {
      this.baseFreqPosThree[0] =
        this.baseFreqPosThree[0] + this.codons[i].fraction;
    } else if (this.codons[i].sequence.charAt(2) == "a") {
      this.baseFreqPosThree[1] =
        this.baseFreqPosThree[1] + this.codons[i].fraction;
    } else if (this.codons[i].sequence.charAt(2) == "t") {
      this.baseFreqPosThree[2] =
        this.baseFreqPosThree[2] + this.codons[i].fraction;
    } else if (this.codons[i].sequence.charAt(2) == "c") {
      this.baseFreqPosThree[3] =
        this.baseFreqPosThree[3] + this.codons[i].fraction;
    }
  }
  return true;
}

//class AminoAcid method fixFraction()
//added to address bug in http://www.kazusa.or.jp/codon/ that causes fraction values to all be given as 0.
function fixFraction() {
  var perThouTotal = 0;
  for (var i = 0; i < this.codons.length; i++) {
    perThouTotal = perThouTotal + this.codons[i].perThou;
  }

  if (perThouTotal == 0) {
    return false;
  }

  for (var i = 0; i < this.codons.length; i++) {
    this.codons[i].fraction = this.codons[i].perThou / perThouTotal;
  }
  return true;
}

//class AminoAcid method getMostCommonCodon()
function setMostCommonCodon() {
  var highestFraction = 0;
  var highestCodon = "nnn";
  for (var i = 0; i < this.codons.length; i++) {
    if (this.codons[i].fraction > highestFraction) {
      highestFraction = this.codons[i].fraction;
      highestCodon = this.codons[i].sequence;
    }
  }
  this.mostCommonCodon = highestCodon;
}

//class AminoAcid method getDegenCodon()
function setDegenCodon() {
  this.degenCodon =
    getConsensusBase(this.baseFreqPosOne) +
    getConsensusBase(this.baseFreqPosTwo) +
    getConsensusBase(this.baseFreqPosThree);
}

//class AminoAcid method getConsensusBase()
function getConsensusBase(baseFreq) {
  var g;
  var a;
  var t;
  var c;
  if (baseFreq[0] > 0) {
    g = true;
  }
  if (baseFreq[1] > 0) {
    a = true;
  }
  if (baseFreq[2] > 0) {
    t = true;
  }
  if (baseFreq[3] > 0) {
    c = true;
  }
  if (!g && !a && !c && !t) {
    return "n";
  }
  if (g && a && c && t) {
    return "n";
  } else if (a && c && t) {
    return "h";
  } else if (a && g && t) {
    return "d";
  } else if (c && g && t) {
    return "b";
  } else if (a && c) {
    return "m";
  } else if (g && t) {
    return "k";
  } else if (a && t) {
    return "w";
  } else if (g && c) {
    return "s";
  } else if (c && t) {
    return "y";
  } else if (a && g) {
    return "r";
  } else if (t) {
    return "t";
  } else if (g) {
    return "g";
  } else if (c) {
    return "c";
  } else if (a) {
    return "a";
  }
  return true;
}

//class AminoAcid
function AminoAcid() {
  this.codons = new Array();
  //g, a, t, c
  this.baseFreqPosOne = new Array(0, 0, 0, 0);
  this.baseFreqPosTwo = new Array(0, 0, 0, 0);
  this.baseFreqPosThree = new Array(0, 0, 0, 0);
  this.mostComonCodon;
  this.degenCodon;
  this.rulerPosOne;
  this.rulerPosTwo;
  this.rulerPosThree;
}

//create and throw away a prototype object
new AminoAcid();

// define object methods
AminoAcid.prototype.fillRuler = fillRuler;
AminoAcid.prototype.getRuler = getRuler;
AminoAcid.prototype.addCodon = addCodon;
AminoAcid.prototype.determineBaseFreq = determineBaseFreq;
AminoAcid.prototype.fixFraction = fixFraction;
AminoAcid.prototype.setMostCommonCodon = setMostCommonCodon;
AminoAcid.prototype.setDegenCodon = setDegenCodon;
AminoAcid.prototype.getConsensusBase = getConsensusBase;

//class Codon
function Codon(sequence, number, perThou, fraction) {
  this.sequence = sequence.toLowerCase();
  this.number = number;
  this.perThou = perThou;
  this.fraction = fraction;
}
