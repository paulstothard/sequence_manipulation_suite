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

function codonPlot(theDocument) {
  var newDna = "";
  var maxInput = 50000000;
  var codonTable;
  var title;

  if (testScript() == false) {
    return false;
  }

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

  newDna = getSequenceFromFasta(theDocument.forms[0].elements[0].value);
  title = getTitleFromFasta(theDocument.forms[0].elements[0].value);
  verifyDna(newDna);
  newDna = removeNonDna(newDna);

  openWindow("Codon Plot");
  outputWindow.document.write(getInfoFromTitleAndSequence(title, newDna));
  openPre();

  writeCodonPlot(codonTable, newDna);
  closePre();
  closeWindow();
  return true;
}

function writeCodonPlot(codonTable, sequence) {
  var markString =
    "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";
  var codon;
  var perThou;
  var fraction;
  var yValue;
  var aminoAcid;

  //replace 'u' with 't'
  sequence = sequence.replace(/u/gi, "t");
  sequence = sequence.replace(/(...)/g, function (str, p1, offset, s) {
    try {
      aminoAcid = codonTable[p1.toString().toLowerCase()].aminoAcid;
      yValue = codonTable[p1.toString().toLowerCase()].fraction;
    } catch (e) {
      aminoAcid = "???";
      yValue = 0;
    }
    return (
      "<b>" +
      p1.toString().toLowerCase() +
      ", " +
      (offset + 1) +
      " to " +
      (offset + 3) +
      " (" +
      aminoAcid +
      ")</b>\n" +
      markString.substring(0, Math.round(yValue * markString.length)) +
      " " +
      yValue.toFixed(2) +
      "\n\n"
    );
  });

  outputWindow.document.write(sequence + "\n");

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
        codonTable[matchArray[2].toLowerCase()].fillCodon(
          matchArray[1],
          parseFloat(matchArray[3]),
          parseFloat(matchArray[4]),
          parseFloat(matchArray[5])
        );
        codonTable.codons.push(matchArray[2].toLowerCase());
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
            matchArray[5]
        );
        return false;
      }
    }
  }

  codonTable.fixFraction();

  return codonTable;
}

//class CodonTable
function CodonTable() {
  this.codons = new Array();
  this.ggg = new Codon();
  this.gga = new Codon();
  this.ggt = new Codon();
  this.ggc = new Codon();
  this.gag = new Codon();
  this.gaa = new Codon();
  this.gat = new Codon();
  this.gac = new Codon();
  this.gtg = new Codon();
  this.gta = new Codon();
  this.gtt = new Codon();
  this.gtc = new Codon();
  this.gcg = new Codon();
  this.gca = new Codon();
  this.gct = new Codon();
  this.gcc = new Codon();
  this.agg = new Codon();
  this.aga = new Codon();
  this.agt = new Codon();
  this.agc = new Codon();
  this.aag = new Codon();
  this.aaa = new Codon();
  this.aat = new Codon();
  this.aac = new Codon();
  this.atg = new Codon();
  this.ata = new Codon();
  this.att = new Codon();
  this.atc = new Codon();
  this.acg = new Codon();
  this.aca = new Codon();
  this.act = new Codon();
  this.acc = new Codon();
  this.tgg = new Codon();
  this.tga = new Codon();
  this.tgt = new Codon();
  this.tgc = new Codon();
  this.tag = new Codon();
  this.taa = new Codon();
  this.tat = new Codon();
  this.tac = new Codon();
  this.ttg = new Codon();
  this.tta = new Codon();
  this.ttt = new Codon();
  this.ttc = new Codon();
  this.tcg = new Codon();
  this.tca = new Codon();
  this.tct = new Codon();
  this.tcc = new Codon();
  this.cgg = new Codon();
  this.cga = new Codon();
  this.cgt = new Codon();
  this.cgc = new Codon();
  this.cag = new Codon();
  this.caa = new Codon();
  this.cat = new Codon();
  this.cac = new Codon();
  this.ctg = new Codon();
  this.cta = new Codon();
  this.ctt = new Codon();
  this.ctc = new Codon();
  this.ccg = new Codon();
  this.cca = new Codon();
  this.cct = new Codon();
  this.ccc = new Codon();
}

//class CodonTable method fixFraction()
//added to address bug in http://www.kazusa.or.jp/codon/ that causes fraction values to all be given as 0.
function fixFraction() {
  for (var i = 0; i < this.codons.length; i++) {
    var outerCodon = this.codons[i];
    var perThouTotal = 0;
    for (var j = 0; j < this.codons.length; j++) {
      var innerCodon = this.codons[j];
      if (this[outerCodon].aminoAcid == this[innerCodon].aminoAcid) {
        perThouTotal = perThouTotal + this[innerCodon].perThou;
      }
    }
    this[outerCodon].fraction = this[outerCodon].perThou / perThouTotal;
  }
  return true;
}

//create and throw away a prototype object
new CodonTable();

// define object methods
CodonTable.prototype.fixFraction = fixFraction;

//class Codon method fillCodon()
function fillCodon(aminoAcid, number, perThou, fraction) {
  this.aminoAcid = aminoAcid;
  this.number = number;
  this.perThou = perThou;
  this.fraction = fraction;
}

//class Codon
function Codon() {
  this.aminoAcid;
  this.number;
  this.perThou;
  this.fraction;
}

//create and throw away a prototype object
new Codon();

// define object methods
Codon.prototype.fillCodon = fillCodon;
