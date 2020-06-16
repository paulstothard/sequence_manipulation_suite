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

function multiRevTrans(theDocument) {
  var newProtein = "";
  var maxInput = 20000000;
  var codonTable;
  var alignArray = new Array();
  var titleArray = new Array();
  var sequenceArray = new Array();

  if (testScript() == false) {
    return false;
  }

  if (
    checkFormElement(theDocument.forms[0].elements[0]) == false ||
    checkTextLength(theDocument.forms[0].elements[0].value, maxInput) ==
      false ||
    checkCodonTable(theDocument.forms[0].elements[4].value) == false
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

  codonTable = makeCodonTable(theDocument.forms[0].elements[4].value);
  if (codonTable == false) {
    return false;
  }
  openWindow("Multi Rev Trans");
  for (var i = 0; i < titleArray.length; i++) {
    outputWindow.document.write(
      getInfoFromTitleAndSequence(titleArray[i], sequenceArray[i])
    );
    if (i < titleArray.length - 1) {
      outputWindow.document.write('<div class="info">Combined with:</div>\n');
    }
  }
  openPre();
  writeConsensusSeq(sequenceArray, codonTable);
  outputWindow.document.write("\n");
  writeMultiRevTrans(sequenceArray, codonTable);
  closePre();
  closeWindow();
  return true;
}

function writeConsensusSeq(sequenceArray, codonTable) {
  var consensusSeq = new Array();
  var aminoAcid;
  var firstG;
  var firstA;
  var firstT;
  var firstC;

  var secondG;
  var secondA;
  var secondT;
  var secondC;

  var thirdG;
  var thirdA;
  var thirdT;
  var thirdC;

  //go through each position in the alignment and fill 12 variables with the sum of the
  for (var i = 0; i < sequenceArray[0].length; i++) {
    firstG = 0;
    firstA = 0;
    firstT = 0;
    firstC = 0;

    secondG = 0;
    secondA = 0;
    secondT = 0;
    secondC = 0;

    thirdG = 0;
    thirdA = 0;
    thirdT = 0;
    thirdC = 0;

    for (var j = 0; j < sequenceArray.length; j++) {
      if (
        sequenceArray[j].charAt(i) == "-" ||
        sequenceArray[j].charAt(i) == "."
      ) {
        firstG = firstG + 0.25;
        firstA = firstA + 0.25;
        firstT = firstT + 0.25;
        firstC = firstC + 0.25;

        secondG = secondG + 0.25;
        secondA = secondA + 0.25;
        secondT = secondT + 0.25;
        secondC = secondC + 0.25;

        thirdG = thirdG + 0.25;
        thirdA = thirdA + 0.25;
        thirdT = thirdT + 0.25;
        thirdC = thirdC + 0.25;
      } else {
        try {
          aminoAcid =
            codonTable[sequenceArray[j].charAt(i).toString().toLowerCase()];
        } catch (e) {
          alert(
            "A codon table entry wasn't found for " + sequenceArray[j].charAt(i)
          );
          return false;
        }

        firstG = firstG + aminoAcid.baseFreqPosOne[0];
        firstA = firstA + aminoAcid.baseFreqPosOne[1];
        firstT = firstT + aminoAcid.baseFreqPosOne[2];
        firstC = firstC + aminoAcid.baseFreqPosOne[3];

        secondG = secondG + aminoAcid.baseFreqPosTwo[0];
        secondA = secondA + aminoAcid.baseFreqPosTwo[1];
        secondT = secondT + aminoAcid.baseFreqPosTwo[2];
        secondC = secondC + aminoAcid.baseFreqPosTwo[3];

        thirdG = thirdG + aminoAcid.baseFreqPosThree[0];
        thirdA = thirdA + aminoAcid.baseFreqPosThree[1];
        thirdT = thirdT + aminoAcid.baseFreqPosThree[2];
        thirdC = thirdC + aminoAcid.baseFreqPosThree[3];
      }
    }

    consensusSeq.push(_getConsensusBase([firstG, firstA, firstT, firstC]));
    consensusSeq.push(_getConsensusBase([secondG, secondA, secondT, secondC]));
    consensusSeq.push(_getConsensusBase([thirdG, thirdA, thirdT, thirdC]));
  }

  outputWindow.document.write(
    "&gt;" +
      "reverse translation of alignment to a " +
      consensusSeq.length +
      " base sequence of consensus codons.\n"
  );
  outputWindow.document.write(addReturns(consensusSeq.join("")));
  outputWindow.document.write("\n");

  return true;
}

function writeMultiRevTrans(sequenceArray, codonTable) {
  var markG =
    "gggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg";
  var markA =
    "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
  var markT =
    "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
  var markC =
    "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";

  var markLength = markG.length;

  var aminoAcid;
  var columnSeq;
  var firstG;
  var firstA;
  var firstT;
  var firstC;

  var secondG;
  var secondA;
  var secondT;
  var secondC;

  var thirdG;
  var thirdA;
  var thirdT;
  var thirdC;

  //go through each position in the alignment and fill 12 variables with the sum of the
  for (var i = 0; i < sequenceArray[0].length; i++) {
    columnSeq = "";

    firstG = 0;
    firstA = 0;
    firstT = 0;
    firstC = 0;

    secondG = 0;
    secondA = 0;
    secondT = 0;
    secondC = 0;

    thirdG = 0;
    thirdA = 0;
    thirdT = 0;
    thirdC = 0;

    for (var j = 0; j < sequenceArray.length; j++) {
      columnSeq = columnSeq + sequenceArray[j].charAt(i);
      if (
        sequenceArray[j].charAt(i) == "-" ||
        sequenceArray[j].charAt(i) == "."
      ) {
        firstG = firstG + 0.25;
        firstA = firstA + 0.25;
        firstT = firstT + 0.25;
        firstC = firstC + 0.25;

        secondG = secondG + 0.25;
        secondA = secondA + 0.25;
        secondT = secondT + 0.25;
        secondC = secondC + 0.25;

        thirdG = thirdG + 0.25;
        thirdA = thirdA + 0.25;
        thirdT = thirdT + 0.25;
        thirdC = thirdC + 0.25;
      } else {
        try {
          aminoAcid =
            codonTable[sequenceArray[j].charAt(i).toString().toLowerCase()];
        } catch (e) {
          alert(
            "A codon table entry wasn't found for " + sequenceArray[j].charAt(i)
          );
          return false;
        }

        firstG = firstG + aminoAcid.baseFreqPosOne[0];
        firstA = firstA + aminoAcid.baseFreqPosOne[1];
        firstT = firstT + aminoAcid.baseFreqPosOne[2];
        firstC = firstC + aminoAcid.baseFreqPosOne[3];

        secondG = secondG + aminoAcid.baseFreqPosTwo[0];
        secondA = secondA + aminoAcid.baseFreqPosTwo[1];
        secondT = secondT + aminoAcid.baseFreqPosTwo[2];
        secondC = secondC + aminoAcid.baseFreqPosTwo[3];

        thirdG = thirdG + aminoAcid.baseFreqPosThree[0];
        thirdA = thirdA + aminoAcid.baseFreqPosThree[1];
        thirdT = thirdT + aminoAcid.baseFreqPosThree[2];
        thirdC = thirdC + aminoAcid.baseFreqPosThree[3];
      }
    }

    firstG = Math.round((markLength * firstG) / sequenceArray.length);
    firstA = Math.round((markLength * firstA) / sequenceArray.length);
    firstT = Math.round((markLength * firstT) / sequenceArray.length);
    firstC = Math.round((markLength * firstC) / sequenceArray.length);

    secondG = Math.round((markLength * secondG) / sequenceArray.length);
    secondA = Math.round((markLength * secondA) / sequenceArray.length);
    secondT = Math.round((markLength * secondT) / sequenceArray.length);
    secondC = Math.round((markLength * secondC) / sequenceArray.length);

    thirdG = Math.round((markLength * thirdG) / sequenceArray.length);
    thirdA = Math.round((markLength * thirdA) / sequenceArray.length);
    thirdT = Math.round((markLength * thirdT) / sequenceArray.length);
    thirdC = Math.round((markLength * thirdC) / sequenceArray.length);

    outputWindow.document.write(
      "<b>" + (i + 1) + "_" + columnSeq + "_" + "first</b>\n"
    );
    outputWindow.document.write(
      "g" +
        markG.substring(0, firstG) +
        " " +
        (firstG / markLength).toFixed(2) +
        "\n" +
        "a" +
        markA.substring(0, firstA) +
        " " +
        (firstA / markLength).toFixed(2) +
        "\n" +
        "T" +
        markT.substring(0, firstT) +
        " " +
        (firstT / markLength).toFixed(2) +
        "\n" +
        "C" +
        markC.substring(0, firstC) +
        " " +
        (firstC / markLength).toFixed(2) +
        "\n"
    );

    outputWindow.document.write(
      "<b>" + (i + 1) + "_" + columnSeq + "_" + "second</b>\n"
    );
    outputWindow.document.write(
      "g" +
        markG.substring(0, secondG) +
        " " +
        (secondG / markLength).toFixed(2) +
        "\n" +
        "a" +
        markA.substring(0, secondA) +
        " " +
        (secondA / markLength).toFixed(2) +
        "\n" +
        "T" +
        markT.substring(0, secondT) +
        " " +
        (secondT / markLength).toFixed(2) +
        "\n" +
        "C" +
        markC.substring(0, secondC) +
        " " +
        (secondC / markLength).toFixed(2) +
        "\n"
    );

    outputWindow.document.write(
      "<b>" + (i + 1) + "_" + columnSeq + "_" + "third</b>\n"
    );
    outputWindow.document.write(
      "g" +
        markG.substring(0, thirdG) +
        " " +
        (thirdG / markLength).toFixed(2) +
        "\n" +
        "a" +
        markA.substring(0, thirdA) +
        " " +
        (thirdA / markLength).toFixed(2) +
        "\n" +
        "T" +
        markT.substring(0, thirdT) +
        " " +
        (thirdT / markLength).toFixed(2) +
        "\n" +
        "C" +
        markC.substring(0, thirdC) +
        " " +
        (thirdC / markLength).toFixed(2) +
        "\n"
    );

    outputWindow.document.write("\n");
  }
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

//class AminoAcid
function AminoAcid() {
  this.codons = new Array();
  this.baseFreqPosOne = new Array(0, 0, 0, 0);
  this.baseFreqPosTwo = new Array(0, 0, 0, 0);
  this.baseFreqPosThree = new Array(0, 0, 0, 0);
  this.rulerPosOne;
  this.rulerPosTwo;
  this.rulerPosThree;
}

//create and throw away a prototype object
new AminoAcid();

// define object methods
AminoAcid.prototype.addCodon = addCodon;
AminoAcid.prototype.determineBaseFreq = determineBaseFreq;
AminoAcid.prototype.fixFraction = fixFraction;

//class Codon
function Codon(sequence, number, perThou, fraction) {
  this.sequence = sequence.toLowerCase();
  this.number = number;
  this.perThou = perThou;
  this.fraction = fraction;
}

function _getConsensusBase(baseFreq) {
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
  return "?";
}
