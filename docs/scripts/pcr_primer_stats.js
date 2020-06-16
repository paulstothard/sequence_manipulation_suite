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

function pcrPrimerStats(theDocument) {
  var newDna = "";
  var title = "";
  var maxInput = 5000000;
  var maxPrimerLength = 50;
  var milliMolarSalt = 50;
  var milliMolarMagnesium = 1.5;
  var nanoMolarPrimerTotal = 200;
  var isPhosphorylated = false;

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

  molarSalt = milliMolarSalt / 1e3; //convert to molar from millimolar
  molarMagnesium = milliMolarMagnesium / 1e3; //convert to molar from millimolar
  molarPrimerTotal = nanoMolarPrimerTotal / 1e9; //convert to molar from nanomolar
  isPhosphorylated = false;

  //molarSalt affects Salt adjusted Tm and Tm (Nearest neighbor)
  //molarMagnesium affects Tm (Nearest neighbor)
  //molarPrimerTotal affects Tm (Nearest neighbor)
  //isPhosphorylated affects molecular weight

  openWindow("PCR Primer Stats");
  openPre();
  var arrayOfFasta = getArrayOfFasta(theDocument.forms[0].elements[0].value);

  outputWindow.document.write("Global settings:\n");
  if (isPhosphorylated) {
    outputWindow.document.write("-The primers have a 5'-phosphate group.\n");
  } else {
    outputWindow.document.write(
      "-The primers do not have a 5'-phosphate group.\n"
    );
  }
  outputWindow.document.write(
    "-Combined concentration of K+ and Na+ in the reaction = " +
      milliMolarSalt +
      " millimolar.\n"
  );
  outputWindow.document.write(
    "-Mg+2 concentration in the reaction = " +
      milliMolarMagnesium +
      " millimolar.\n"
  );
  outputWindow.document.write(
    "-Primer concentration in the reaction = " +
      nanoMolarPrimerTotal +
      " nanomolar.\n"
  );
  outputWindow.document.write("\n");
  outputWindow.document.write(
    "------------------------------------------------------------\n"
  );

  for (var i = 0; i < arrayOfFasta.length; i++) {
    newDna = getSequenceFromFasta(arrayOfFasta[i]);
    title = getTitleFromFasta(arrayOfFasta[i]);

    newDna = _removeNonPrimer(newDna);

    if (newDna.length == 0) {
      continue;
    }

    if (newDna.length > maxPrimerLength) {
      continue;
    }

    var percentGC = _percentGC(newDna);
    var nearestNeighborTm = _nearestNeighborTm(
      newDna,
      molarSalt,
      molarPrimerTotal,
      molarMagnesium
    );
    var selfCompHash = _getSelfComplementarityReport(newDna, 3, 50);
    var hairpinHash = _getHairpinReport(newDna, 3, 50);

    outputWindow.document.write(
      "------------------------------------------------------------\n"
    );
    outputWindow.document.write("General properties:\n");
    outputWindow.document.write("-------------------\n");
    outputWindow.document.write(
      rightNum("Primer name:", "", 32, "") + title + "\n"
    );
    outputWindow.document.write(
      rightNum("Primer sequence:", "", 32, "") + newDna + "\n"
    );
    outputWindow.document.write(
      rightNum("Sequence length:", "", 32, "") + newDna.length + "\n"
    );
    outputWindow.document.write(
      rightNum("Base counts:", "", 32, "") + _baseCounts(newDna) + "\n"
    );
    outputWindow.document.write(
      rightNum("GC content (%):", "", 32, "") + percentGC + "\n"
    );
    outputWindow.document.write(
      rightNum("Molecular weight (Daltons):", "", 32, "") +
        _molecularWeight(newDna, isPhosphorylated) +
        "\n"
    );
    outputWindow.document.write(
      rightNum("nmol/A260:", "", 32, "") + _nmolPerA260(newDna) + "\n"
    );
    outputWindow.document.write(
      rightNum("micrograms/A260:", "", 32, "") +
        _microgramsPerA260(newDna, isPhosphorylated) +
        "\n"
    );
    outputWindow.document.write(
      rightNum("Basic Tm (degrees C):", "", 32, "") + _basicTm(newDna) + "\n"
    );
    outputWindow.document.write(
      rightNum("Salt adjusted Tm (degrees C):", "", 32, "") +
        _molarSaltAdjustedTm(newDna, molarSalt) +
        "\n"
    );
    outputWindow.document.write(
      rightNum("Nearest neighbor Tm (degrees C):", "", 32, "") +
        nearestNeighborTm +
        "\n"
    );

    outputWindow.document.write("\n");
    outputWindow.document.write("PCR suitability tests (Pass / Warning):\n");
    outputWindow.document.write("------------------------------------\n");
    outputWindow.document.write(
      rightNum("Single base runs:", "", 32, "") +
        _getBaseRunsReport(newDna, 5) +
        "\n"
    );
    outputWindow.document.write(
      rightNum("Dinucleotide base runs:", "", 32, "") +
        _getDiNucleotideRunsReport(newDna, 5) +
        "\n"
    );
    outputWindow.document.write(
      rightNum("Length:", "", 32, "") +
        _getSuitableLengthReport(newDna, 14, 30) +
        "\n"
    );
    outputWindow.document.write(
      rightNum("Percent GC:", "", 32, "") +
        _getSuitableGCReport(newDna, percentGC, 40, 60) +
        "\n"
    );
    outputWindow.document.write(
      rightNum("Tm (Nearest neighbor):", "", 32, "") +
        _getSuitableTmReport(newDna, nearestNeighborTm, 50, 58) +
        "\n"
    );
    outputWindow.document.write(
      rightNum("GC clamp:", "", 32, "") +
        _getSuitableThreePrimeGC(newDna, 1, 3) +
        "\n"
    );
    outputWindow.document.write(
      rightNum("Self-annealing:", "", 32, "") + selfCompHash["report"] + "\n"
    );
    if (selfCompHash["report"] != "Pass") {
      outputWindow.document.write(
        rightNum(":", "", 32, "") + selfCompHash["upper"] + "\n"
      );
      outputWindow.document.write(
        rightNum(":", "", 32, "") + selfCompHash["divider"] + "\n"
      );
      outputWindow.document.write(
        rightNum(":", "", 32, "") + selfCompHash["lower"] + "\n"
      );
    }
    outputWindow.document.write(
      rightNum("Hairpin formation:", "", 32, "") + hairpinHash["report"] + "\n"
    );
    if (hairpinHash["report"] != "Pass") {
      outputWindow.document.write(
        rightNum(":", "", 32, "") + hairpinHash["upper"] + "\n"
      );
      outputWindow.document.write(
        rightNum(":", "", 32, "") + hairpinHash["divider"] + "\n"
      );
      outputWindow.document.write(
        rightNum(":", "", 32, "") + hairpinHash["lower"] + "\n"
      );
    }
    outputWindow.document.write(
      "------------------------------------------------------------\n"
    );
    outputWindow.document.write("\n");
  }
  closePre();
  closeWindow();
  return true;
}

//Nucleotide Code: Base:
//---------------- -----
//A.................Adenine
//C.................Cytosine
//G.................Guanine
//T (or U)..........Thymine (or Uracil)
//R.................A or G
//Y.................C or T
//S.................G or C
//W.................A or T
//K.................G or T
//M.................A or C
//B.................C or G or T
//D.................A or G or T
//H.................A or C or T
//V.................A or C or G
//N.................any base

function _removeNonPrimer(sequence) {
  sequence.replace(/u/g, "t");
  sequence.replace(/U/g, "T");
  return sequence.replace(/[^gatcryswkmbdhvnGATCRYSWKMBDHVN]/g, "");
}

function _containsOnlyNonDegenerates(sequence) {
  if (sequence.search(/[^gatc]/i) == -1) {
    return true;
  }
  return false;
}

function _baseCounts(sequence) {
  var numG = _getBaseCount(sequence, "g");
  var numA = _getBaseCount(sequence, "a");
  var numT = _getBaseCount(sequence, "t");
  var numC = _getBaseCount(sequence, "c");
  var numOther = sequence.length - (numG + numA + numT + numC);
  return (
    "G=" +
    numG +
    "; A=" +
    numA +
    "; T=" +
    numT +
    "; C=" +
    numC +
    "; Other=" +
    numOther +
    ";"
  );
}

function _microgramsPerA260(sequence, isPhosphorylated) {
  if (_containsOnlyNonDegenerates(sequence)) {
    return _microgramsPerA260NonDegen(sequence, isPhosphorylated);
  } else {
    return _microgramsPerA260Degen(sequence, isPhosphorylated);
  }
}

function _microgramsPerA260NonDegen(sequence, isPhosphorylated) {
  var mw = _mw(sequence, isPhosphorylated);
  var result = mw / _getExtinctionCoefficient(sequence);
  return result.toFixed(2);
}

function _microgramsPerA260Degen(sequence, isPhosphorylated) {
  var lowerBoundsSequence = sequence;
  var upperBoundsSequence = sequence;

  //replace all other degenerates with the base with lowest value in lowerBoundsSequence
  lowerBoundsSequence = lowerBoundsSequence.replace(/r/gi, "g");
  lowerBoundsSequence = lowerBoundsSequence.replace(/y/gi, "c");
  lowerBoundsSequence = lowerBoundsSequence.replace(/s/gi, "c");
  lowerBoundsSequence = lowerBoundsSequence.replace(/w/gi, "t");
  lowerBoundsSequence = lowerBoundsSequence.replace(/k/gi, "t");
  lowerBoundsSequence = lowerBoundsSequence.replace(/m/gi, "c");
  lowerBoundsSequence = lowerBoundsSequence.replace(/b/gi, "c");
  lowerBoundsSequence = lowerBoundsSequence.replace(/d/gi, "t");
  lowerBoundsSequence = lowerBoundsSequence.replace(/h/gi, "c");
  lowerBoundsSequence = lowerBoundsSequence.replace(/v/gi, "c");
  lowerBoundsSequence = lowerBoundsSequence.replace(/n/gi, "c");

  //replace all other degenerates with base with highest value in upperBoundsSequence
  upperBoundsSequence = upperBoundsSequence.replace(/r/gi, "a");
  upperBoundsSequence = upperBoundsSequence.replace(/y/gi, "t");
  upperBoundsSequence = upperBoundsSequence.replace(/s/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/w/gi, "a");
  upperBoundsSequence = upperBoundsSequence.replace(/k/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/m/gi, "a");
  upperBoundsSequence = upperBoundsSequence.replace(/b/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/d/gi, "a");
  upperBoundsSequence = upperBoundsSequence.replace(/h/gi, "a");
  upperBoundsSequence = upperBoundsSequence.replace(/v/gi, "a");
  upperBoundsSequence = upperBoundsSequence.replace(/n/gi, "a");

  //swap upper and lower because of how downstream calculation is done
  //return _microgramsPerA260NonDegen(lowerBoundsSequence, isPhosphorylated) + " to " + _microgramsPerA260NonDegen(upperBoundsSequence, isPhosphorylated);

  return (
    _microgramsPerA260NonDegen(upperBoundsSequence, isPhosphorylated) +
    " to " +
    _microgramsPerA260NonDegen(lowerBoundsSequence, isPhosphorylated)
  );
}

function _nmolPerA260(sequence) {
  if (_containsOnlyNonDegenerates(sequence)) {
    return _nmolPerA260NonDegen(sequence);
  } else {
    return _nmolPerA260Degen(sequence);
  }
}

function _nmolPerA260NonDegen(sequence) {
  var result = _getExtinctionCoefficient(sequence);
  result = (1 / result) * 1000;
  return result.toFixed(2);
}

function _getExtinctionCoefficient(sequence) {
  var dimerValues = _getDimerExtinctionCoefficients();
  var singleValues = _getSingleExtinctionCoefficients();
  var dimerSum = 0;
  var singleSum = 0;
  sequence = sequence.toLowerCase();

  for (var i = 0; i < sequence.length - 1; i++) {
    dimerSum =
      dimerSum + dimerValues[sequence.charAt(i) + sequence.charAt(i + 1)];
  }

  for (var i = 1; i < sequence.length - 1; i++) {
    singleSum = singleSum + singleValues[sequence.charAt(i)];
  }
  return 2 * dimerSum - singleSum;
}

function _nmolPerA260Degen(sequence) {
  var lowerBoundsSequence = sequence;
  var upperBoundsSequence = sequence;

  //replace all other degenerates with the base with lowest value in lowerBoundsSequence
  lowerBoundsSequence = lowerBoundsSequence.replace(/r/gi, "g");
  lowerBoundsSequence = lowerBoundsSequence.replace(/y/gi, "c");
  lowerBoundsSequence = lowerBoundsSequence.replace(/s/gi, "c");
  lowerBoundsSequence = lowerBoundsSequence.replace(/w/gi, "t");
  lowerBoundsSequence = lowerBoundsSequence.replace(/k/gi, "t");
  lowerBoundsSequence = lowerBoundsSequence.replace(/m/gi, "c");
  lowerBoundsSequence = lowerBoundsSequence.replace(/b/gi, "c");
  lowerBoundsSequence = lowerBoundsSequence.replace(/d/gi, "t");
  lowerBoundsSequence = lowerBoundsSequence.replace(/h/gi, "c");
  lowerBoundsSequence = lowerBoundsSequence.replace(/v/gi, "c");
  lowerBoundsSequence = lowerBoundsSequence.replace(/n/gi, "c");

  //replace all other degenerates with base with highest value in upperBoundsSequence
  upperBoundsSequence = upperBoundsSequence.replace(/r/gi, "a");
  upperBoundsSequence = upperBoundsSequence.replace(/y/gi, "t");
  upperBoundsSequence = upperBoundsSequence.replace(/s/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/w/gi, "a");
  upperBoundsSequence = upperBoundsSequence.replace(/k/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/m/gi, "a");
  upperBoundsSequence = upperBoundsSequence.replace(/b/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/d/gi, "a");
  upperBoundsSequence = upperBoundsSequence.replace(/h/gi, "a");
  upperBoundsSequence = upperBoundsSequence.replace(/v/gi, "a");
  upperBoundsSequence = upperBoundsSequence.replace(/n/gi, "a");

  //swap upper and lower because of how downstream calculation is done
  //return _nmolPerA260NonDegen(lowerBoundsSequence) + " to " + _nmolPerA260NonDegen(upperBoundsSequence);
  return (
    _nmolPerA260NonDegen(upperBoundsSequence) +
    " to " +
    _nmolPerA260NonDegen(lowerBoundsSequence)
  );
}

function _percentGC(sequence) {
  if (_containsOnlyNonDegenerates(sequence)) {
    return _percentGCNonDegen(sequence);
  } else {
    return _percentGCDegen(sequence);
  }
}

function _percentGCNonDegen(sequence) {
  var numHits = _getBaseCount(sequence, "g") + _getBaseCount(sequence, "c");
  return ((numHits / sequence.length) * 100).toFixed(2);
}

function _percentGCDegen(sequence) {
  var lowerBoundsSequence = sequence;
  var upperBoundsSequence = sequence;

  //replace degenerates that must be g or c with g in both sequences
  lowerBoundsSequence = lowerBoundsSequence.replace(/s/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/s/gi, "g");

  //replace degenerates that must be a or t with a in both sequences
  lowerBoundsSequence = lowerBoundsSequence.replace(/w/gi, "a");
  upperBoundsSequence = upperBoundsSequence.replace(/w/gi, "a");

  //replace all other degenerates with a or t in lowerBoundsSequence
  lowerBoundsSequence = lowerBoundsSequence.replace(/r/gi, "a");
  lowerBoundsSequence = lowerBoundsSequence.replace(/y/gi, "t");
  lowerBoundsSequence = lowerBoundsSequence.replace(/k/gi, "t");
  lowerBoundsSequence = lowerBoundsSequence.replace(/m/gi, "a");
  lowerBoundsSequence = lowerBoundsSequence.replace(/b/gi, "t");
  lowerBoundsSequence = lowerBoundsSequence.replace(/d/gi, "a");
  lowerBoundsSequence = lowerBoundsSequence.replace(/h/gi, "a");
  lowerBoundsSequence = lowerBoundsSequence.replace(/v/gi, "a");
  lowerBoundsSequence = lowerBoundsSequence.replace(/n/gi, "a");

  //replace all other degenerates with g or c in upperBoundsSequence
  upperBoundsSequence = upperBoundsSequence.replace(/r/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/y/gi, "c");
  upperBoundsSequence = upperBoundsSequence.replace(/k/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/m/gi, "c");
  upperBoundsSequence = upperBoundsSequence.replace(/b/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/d/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/h/gi, "c");
  upperBoundsSequence = upperBoundsSequence.replace(/v/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/n/gi, "g");

  return (
    _percentGCNonDegen(lowerBoundsSequence) +
    " to " +
    _percentGCNonDegen(upperBoundsSequence)
  );
}

function _molecularWeight(sequence, isPhosphorylated) {
  if (_containsOnlyNonDegenerates(sequence)) {
    return _molecularWeightNonDegen(sequence, isPhosphorylated);
  } else {
    return _molecularWeightDegen(sequence, isPhosphorylated);
  }
}

function _molecularWeightNonDegen(sequence, isPhosphorylated) {
  return _mw(sequence, isPhosphorylated).toFixed(2);
}

function _mw(sequence, isPhosphorylated) {
  //DNA molecular weight for synthesized oligonucleotides
  var g = _getBaseCount(sequence, "g");
  var a = _getBaseCount(sequence, "a");
  var t = _getBaseCount(sequence, "t");
  var c = _getBaseCount(sequence, "c");
  var phosAdjust = 0;
  if (isPhosphorylated) {
    phosAdjust = 79.0;
  }
  return g * 329.21 + a * 313.21 + t * 304.2 + c * 289.18 - 61.96 + phosAdjust;
}

function _molecularWeightDegen(sequence, isPhosphorylated) {
  var lowerBoundsSequence = sequence;
  var upperBoundsSequence = sequence;

  //replace all other degenerates with lightest base possible in lowerBoundsSequence
  lowerBoundsSequence = lowerBoundsSequence.replace(/r/gi, "a");
  lowerBoundsSequence = lowerBoundsSequence.replace(/y/gi, "c");
  lowerBoundsSequence = lowerBoundsSequence.replace(/s/gi, "c");
  lowerBoundsSequence = lowerBoundsSequence.replace(/w/gi, "t");
  lowerBoundsSequence = lowerBoundsSequence.replace(/k/gi, "t");
  lowerBoundsSequence = lowerBoundsSequence.replace(/m/gi, "c");
  lowerBoundsSequence = lowerBoundsSequence.replace(/b/gi, "c");
  lowerBoundsSequence = lowerBoundsSequence.replace(/d/gi, "t");
  lowerBoundsSequence = lowerBoundsSequence.replace(/h/gi, "c");
  lowerBoundsSequence = lowerBoundsSequence.replace(/v/gi, "c");
  lowerBoundsSequence = lowerBoundsSequence.replace(/n/gi, "c");

  //replace all other degenerates with heaviest base possible in upperBoundsSequence
  upperBoundsSequence = upperBoundsSequence.replace(/r/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/y/gi, "t");
  upperBoundsSequence = upperBoundsSequence.replace(/s/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/w/gi, "a");
  upperBoundsSequence = upperBoundsSequence.replace(/k/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/m/gi, "a");
  upperBoundsSequence = upperBoundsSequence.replace(/b/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/d/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/h/gi, "a");
  upperBoundsSequence = upperBoundsSequence.replace(/v/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/n/gi, "g");

  return (
    _molecularWeightNonDegen(lowerBoundsSequence, isPhosphorylated) +
    " to " +
    _molecularWeightNonDegen(upperBoundsSequence, isPhosphorylated)
  );
}

function _basicTm(sequence) {
  if (_containsOnlyNonDegenerates(sequence)) {
    return _basicTmNonDegen(sequence);
  } else {
    return _basicTmDegen(sequence);
  }
}

function _basicTmNonDegen(sequence) {
  //Simple formula when primer length < 14 bases from
  //Rychlik, W. and Rhoads, R.E. (1989) Nucleic Acids Research 17, 8543
  //Tm = 4C x (number of G's and C's in the primer) + 2C x (number of A's and T's in the primer)
  //
  //When longer use:
  //tm = 64.9C + 41C * (number of G's and C's in the primer - 16.4)/ primer length
  //
  //both assume reaction is carried out in the presence of 50mM monovalent cations

  if (sequence.length < 14) {
    var numG = _getBaseCount(sequence, "g");
    var numC = _getBaseCount(sequence, "c");
    var numA = _getBaseCount(sequence, "a");
    var numT = _getBaseCount(sequence, "t");
    return (4 * (numG + numC) + 2 * (numA + numT)).toFixed(0);
  } else {
    var numG = _getBaseCount(sequence, "g");
    var numC = _getBaseCount(sequence, "c");
    return (64.9 + (41 * (numG + numC - 16.4)) / sequence.length).toFixed(0);
  }
}

function _basicTmDegen(sequence) {
  var lowerBoundsSequence = sequence;
  var upperBoundsSequence = sequence;

  //replace degenerates that must be g or c with g in both sequences
  lowerBoundsSequence = lowerBoundsSequence.replace(/s/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/s/gi, "g");

  //replace degenerates that must be a or t with a in both sequences
  lowerBoundsSequence = lowerBoundsSequence.replace(/w/gi, "a");
  upperBoundsSequence = upperBoundsSequence.replace(/w/gi, "a");

  //replace all other degenerates with a or t in lowerBoundsSequence
  lowerBoundsSequence = lowerBoundsSequence.replace(/r/gi, "a");
  lowerBoundsSequence = lowerBoundsSequence.replace(/y/gi, "t");
  lowerBoundsSequence = lowerBoundsSequence.replace(/k/gi, "t");
  lowerBoundsSequence = lowerBoundsSequence.replace(/m/gi, "a");
  lowerBoundsSequence = lowerBoundsSequence.replace(/b/gi, "t");
  lowerBoundsSequence = lowerBoundsSequence.replace(/d/gi, "a");
  lowerBoundsSequence = lowerBoundsSequence.replace(/h/gi, "a");
  lowerBoundsSequence = lowerBoundsSequence.replace(/v/gi, "a");
  lowerBoundsSequence = lowerBoundsSequence.replace(/n/gi, "a");

  //replace all other degenerates with g or c in upperBoundsSequence
  upperBoundsSequence = upperBoundsSequence.replace(/r/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/y/gi, "c");
  upperBoundsSequence = upperBoundsSequence.replace(/k/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/m/gi, "c");
  upperBoundsSequence = upperBoundsSequence.replace(/b/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/d/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/h/gi, "c");
  upperBoundsSequence = upperBoundsSequence.replace(/v/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/n/gi, "g");

  return (
    _basicTmNonDegen(lowerBoundsSequence) +
    " to " +
    _basicTmNonDegen(upperBoundsSequence)
  );
}

//molarSalt in molar concentration
function _molarSaltAdjustedTm(sequence, molarSalt) {
  if (_containsOnlyNonDegenerates(sequence)) {
    return _molarSaltAdjustedTmNonDegen(sequence, molarSalt);
  } else {
    return _molarSaltAdjustedTmDegen(sequence, molarSalt);
  }
}

function _molarSaltAdjustedTmNonDegen(sequence, molarSalt) {
  //commonly used formula takes into account the molarSalt concentration of the reaction:
  //Tm = 81.5C + 7.21C x Math.log(molarSalt) + (0.41 x GC) - (675 / primer length);
  //see refs
  //Rychlik, W. and Rhoads, R.E. (1989) Nucl. Acids Res. 17, 8543.
  //PCR Core Systems Technical Bulletin #TB254, Promega Corporation.
  //Sambrook, J., Fritsch, E.F. and Maniatis, T. (1989) Molecular Cloning: A Laboratory Manual, Cold Spring Harbor Laboratory Press, Cold Spring Harbor, NY.
  //Mueller, P.R. et al. (1993) In: Current Protocols in Molecular Biology 15.5, Greene Publishing Associates, Inc. and John Wiley and Sons, New York.

  var gcHits = _getBaseCount(sequence, "g") + _getBaseCount(sequence, "c");
  var pGC = (gcHits / sequence.length) * 100;
  return (
    81.5 +
    7.21 * Math.log(molarSalt) +
    0.41 * pGC -
    675 / sequence.length
  ).toFixed(0);
}

function _molarSaltAdjustedTmDegen(sequence, molarSalt) {
  var lowerBoundsSequence = sequence;
  var upperBoundsSequence = sequence;

  //replace degenerates that must be g or c with g in both sequences
  lowerBoundsSequence = lowerBoundsSequence.replace(/s/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/s/gi, "g");

  //replace degenerates that must be a or t with a in both sequences
  lowerBoundsSequence = lowerBoundsSequence.replace(/w/gi, "a");
  upperBoundsSequence = upperBoundsSequence.replace(/w/gi, "a");

  //replace all other degenerates with a or t in lowerBoundsSequence
  lowerBoundsSequence = lowerBoundsSequence.replace(/r/gi, "a");
  lowerBoundsSequence = lowerBoundsSequence.replace(/y/gi, "t");
  lowerBoundsSequence = lowerBoundsSequence.replace(/k/gi, "t");
  lowerBoundsSequence = lowerBoundsSequence.replace(/m/gi, "a");
  lowerBoundsSequence = lowerBoundsSequence.replace(/b/gi, "t");
  lowerBoundsSequence = lowerBoundsSequence.replace(/d/gi, "a");
  lowerBoundsSequence = lowerBoundsSequence.replace(/h/gi, "a");
  lowerBoundsSequence = lowerBoundsSequence.replace(/v/gi, "a");
  lowerBoundsSequence = lowerBoundsSequence.replace(/n/gi, "a");

  //replace all other degenerates with g or c in upperBoundsSequence
  upperBoundsSequence = upperBoundsSequence.replace(/r/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/y/gi, "c");
  upperBoundsSequence = upperBoundsSequence.replace(/k/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/m/gi, "c");
  upperBoundsSequence = upperBoundsSequence.replace(/b/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/d/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/h/gi, "c");
  upperBoundsSequence = upperBoundsSequence.replace(/v/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/n/gi, "g");

  return (
    _molarSaltAdjustedTmNonDegen(lowerBoundsSequence, molarSalt) +
    " to " +
    _molarSaltAdjustedTmNonDegen(upperBoundsSequence, molarSalt)
  );
}

function _nearestNeighborTm(
  sequence,
  molarSalt,
  molarPrimerTotal,
  molarMagnesium
) {
  if (_containsOnlyNonDegenerates(sequence)) {
    return _nearestNeighborTmNonDegen(
      sequence,
      molarSalt,
      molarPrimerTotal,
      molarMagnesium
    );
  } else {
    return _nearestNeighborTmDegen(
      sequence,
      molarSalt,
      molarPrimerTotal,
      molarMagnesium
    );
  }
}

function _nearestNeighborTmNonDegen(
  sequence,
  molarSalt,
  molarPrimerTotal,
  molarMagnesium
) {
  //The most sophisticated Tm calculations take into account the exact sequence and base stacking parameters, not just the base composition.
  //Tm = ((1000* dh)/(ds+(R * Math.log(primer concentration))))-273.15;
  //Borer P.N. et al. (1974)  J. Mol. Biol. 86, 843.
  //SantaLucia, J. (1998) Proc. Nat. Acad. Sci. USA 95, 1460.
  //Allawi, H.T. and SantaLucia, J. Jr. (1997) Biochemistry 36, 10581.
  //von Ahsen N. et al. (1999) Clin. Chem. 45, 2094.

  sequence = sequence.toLowerCase();

  var R = 1.987; //universal gas constant in Cal/degrees C * mol
  var ds = 0; //cal/Kelvin/mol
  var dh = 0; //kcal/mol

  //perform salt correction
  var correctedSalt = molarSalt + molarMagnesium * 140; //adjust for greater stabilizing effects of Mg compared to Na or K. See von Ahsen et al 1999
  ds = ds + 0.368 * (sequence.length - 1) * Math.log(correctedSalt); //from von Ahsen et al 1999

  //perform terminal corrections
  var termDsCorr = _getTerminalCorrectionsDsHash();
  ds = ds + termDsCorr[sequence.charAt(0)];
  ds = ds + termDsCorr[sequence.charAt(sequence.length - 1)];

  var termDhCorr = _getTerminalCorrectionsDhHash();
  dh = dh + termDhCorr[sequence.charAt(0)];
  dh = dh + termDhCorr[sequence.charAt(sequence.length - 1)];

  var dsValues = _getDsHash();
  var dhValues = _getDhHash();

  for (var i = 0; i < sequence.length - 1; i++) {
    ds = ds + dsValues[sequence.charAt(i) + sequence.charAt(i + 1)];
    dh = dh + dhValues[sequence.charAt(i) + sequence.charAt(i + 1)];
  }
  return (
    (1000 * dh) / (ds + R * Math.log(molarPrimerTotal / 2)) -
    273.15
  ).toFixed(2);
}

function _nearestNeighborTmDegen(
  sequence,
  molarSalt,
  molarPrimerTotal,
  molarMagnesium
) {
  var lowerBoundsSequence = sequence;
  var upperBoundsSequence = sequence;

  //replace degenerates that must be g or c with g in both sequences
  lowerBoundsSequence = lowerBoundsSequence.replace(/s/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/s/gi, "g");

  //replace degenerates that must be a or t with a in both sequences
  lowerBoundsSequence = lowerBoundsSequence.replace(/w/gi, "a");
  upperBoundsSequence = upperBoundsSequence.replace(/w/gi, "a");

  //replace all other degenerates with a or t in lowerBoundsSequence
  lowerBoundsSequence = lowerBoundsSequence.replace(/r/gi, "a");
  lowerBoundsSequence = lowerBoundsSequence.replace(/y/gi, "t");
  lowerBoundsSequence = lowerBoundsSequence.replace(/k/gi, "t");
  lowerBoundsSequence = lowerBoundsSequence.replace(/m/gi, "a");
  lowerBoundsSequence = lowerBoundsSequence.replace(/b/gi, "t");
  lowerBoundsSequence = lowerBoundsSequence.replace(/d/gi, "a");
  lowerBoundsSequence = lowerBoundsSequence.replace(/h/gi, "a");
  lowerBoundsSequence = lowerBoundsSequence.replace(/v/gi, "a");
  lowerBoundsSequence = lowerBoundsSequence.replace(/n/gi, "a");

  //replace all other degenerates with g or c in upperBoundsSequence
  upperBoundsSequence = upperBoundsSequence.replace(/r/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/y/gi, "c");
  upperBoundsSequence = upperBoundsSequence.replace(/k/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/m/gi, "c");
  upperBoundsSequence = upperBoundsSequence.replace(/b/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/d/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/h/gi, "c");
  upperBoundsSequence = upperBoundsSequence.replace(/v/gi, "g");
  upperBoundsSequence = upperBoundsSequence.replace(/n/gi, "g");

  return (
    _nearestNeighborTmNonDegen(
      lowerBoundsSequence,
      molarSalt,
      molarPrimerTotal,
      molarMagnesium
    ) +
    " to " +
    _nearestNeighborTmNonDegen(
      upperBoundsSequence,
      molarSalt,
      molarPrimerTotal,
      molarMagnesium
    )
  );
}

function _getBaseCount(sequence, base) {
  var basePattern = new RegExp(base, "gi");
  if (sequence.search(basePattern) != -1) {
    return sequence.match(basePattern).length;
  } else {
    return 0;
  }
}

function _getTerminalCorrectionsDsHash() {
  //SantaLucia, J. (1998) Proc. Nat. Acad. Sci. USA 95, 1460.
  var hash = {};
  hash["g"] = -2.8;
  hash["a"] = 4.1;
  hash["t"] = 4.1;
  hash["c"] = -2.8;
  return hash;
}

function _getTerminalCorrectionsDhHash() {
  //SantaLucia, J. (1998) Proc. Nat. Acad. Sci. USA 95, 1460.
  var hash = {};
  hash["g"] = 0.1;
  hash["a"] = 2.3;
  hash["t"] = 2.3;
  hash["c"] = 0.1;
  return hash;
}

function _getDsHash() {
  //SantaLucia, J. (1998) Proc. Nat. Acad. Sci. USA 95, 1460.
  var hash = {};
  hash["gg"] = -19.9;
  hash["ga"] = -22.2;
  hash["gt"] = -22.4;
  hash["gc"] = -27.2;

  hash["ag"] = -21.0;
  hash["aa"] = -22.2;
  hash["at"] = -20.4;
  hash["ac"] = -22.4;

  hash["tg"] = -22.7;
  hash["ta"] = -21.3;
  hash["tt"] = -22.2;
  hash["tc"] = -22.2;

  hash["cg"] = -27.2;
  hash["ca"] = -22.7;
  hash["ct"] = -21.0;
  hash["cc"] = -19.9;

  return hash;
}

function _getDhHash() {
  //SantaLucia, J. (1998) Proc. Nat. Acad. Sci. USA 95, 1460.
  var hash = {};
  hash["gg"] = -8.0;
  hash["ga"] = -8.2;
  hash["gt"] = -8.4;
  hash["gc"] = -10.6;

  hash["ag"] = -7.8;
  hash["aa"] = -7.9;
  hash["at"] = -7.2;
  hash["ac"] = -8.4;

  hash["tg"] = -8.5;
  hash["ta"] = -7.2;
  hash["tt"] = -7.9;
  hash["tc"] = -8.2;

  hash["cg"] = -10.6;
  hash["ca"] = -8.5;
  hash["ct"] = -7.8;
  hash["cc"] = -8.0;

  return hash;
}

function _getDimerExtinctionCoefficients() {
  //netprimer documentation
  var hash = {};
  hash["gg"] = 10.8;
  hash["ga"] = 12.6;
  hash["gt"] = 10.0;
  hash["gc"] = 8.8;

  hash["ag"] = 12.5;
  hash["aa"] = 13.7;
  hash["at"] = 11.4;
  hash["ac"] = 10.6;

  hash["tg"] = 9.5;
  hash["ta"] = 11.7;
  hash["tt"] = 8.4;
  hash["tc"] = 8.1;

  hash["cg"] = 9.0;
  hash["ca"] = 10.6;
  hash["ct"] = 7.6;
  hash["cc"] = 7.3;

  return hash;
}

function _getSingleExtinctionCoefficients() {
  //netprimer documentation
  var hash = {};
  hash["g"] = 11.5;
  hash["a"] = 15.4;
  hash["t"] = 8.7;
  hash["c"] = 7.4;

  return hash;
}

function _getBaseRunsReport(sequence, minRunLength) {
  var report = "";
  var hasRun = false;
  var nucleotides = ["G", "A", "T", "C"];

  for (var i = 0; i < nucleotides.length; i++) {
    if (_hasRunOfBases(sequence, nucleotides[i], minRunLength)) {
      hasRun = true;
      report = report + "Contains run of " + nucleotides[i] + "'s; ";
    }
  }

  if (hasRun) {
    return "Warning: " + report;
  } else {
    return "Pass";
  }
}

function _getDiNucleotideRunsReport(sequence, minRunLength) {
  var report = "";
  var hasRun = false;
  var diNucleotides = [
    "GA",
    "GT",
    "GC",
    "AG",
    "AT",
    "AC",
    "TG",
    "TA",
    "TC",
    "CG",
    "CA",
    "CT",
  ];

  for (var i = 0; i < diNucleotides.length; i++) {
    if (_hasRunOfBases(sequence, diNucleotides[i], minRunLength)) {
      hasRun = true;
      report = report + "Contains run of " + diNucleotides[i] + "'s; ";
    }
  }

  if (hasRun) {
    return "Warning: " + report;
  } else {
    return "Pass";
  }
}

function _hasRunOfBases(sequence, base, minRunLength) {
  var basePattern = new RegExp("(?:" + base + "){" + minRunLength + ",}", "gi");
  if (sequence.search(basePattern) != -1) {
    return sequence.match(basePattern).length;
  } else {
    return 0;
  }
}

function _getSuitableLengthReport(
  sequence,
  minSuitableLength,
  maxSuitableLength
) {
  var report = "";
  var hasProblem = false;

  if (sequence.length < minSuitableLength) {
    hasProblem = true;
    report = report + "Contains fewer than " + minSuitableLength + " bases; ";
  }

  if (sequence.length > maxSuitableLength) {
    hasProblem = true;
    report = report + "Contains more than " + maxSuitableLength + " bases; ";
  }

  if (hasProblem) {
    return "Warning: " + report;
  } else {
    return "Pass";
  }
}

function _getSuitableGCReport(
  sequence,
  percentGCRange,
  minSuitableGC,
  maxSuitableGC
) {
  var report = "";
  var hasProblem = false;
  var lowerCalculated;
  var upperCalculated;

  //percentGCRange may be a single number or a string containing something like "40 to 60";
  var rangePattern = new RegExp("([d.]+)D+([d.]+)", "gi");
  if (percentGCRange.search(rangePattern) != -1) {
    lowerCalculated = parseFloat($1);
    upperCalculated = parseFloat($2);
  } else {
    lowerCalculated = parseFloat(percentGCRange);
    upperCalculated = parseFloat(percentGCRange);
  }

  if (lowerCalculated < minSuitableGC) {
    hasProblem = true;
    report = report + "%GC is less than " + minSuitableGC + "; ";
  }

  if (upperCalculated > maxSuitableGC) {
    hasProblem = true;
    report = report + "%GC is greater than " + maxSuitableGC + "; ";
  }

  if (hasProblem) {
    return "Warning: " + report;
  } else {
    return "Pass";
  }
}

function _getSuitableLengthReport(
  sequence,
  minSuitableLength,
  maxSuitableLength
) {
  var report = "";
  var hasProblem = false;

  if (sequence.length < minSuitableLength) {
    hasProblem = true;
    report = report + "Contains fewer than " + minSuitableLength + " bases; ";
  }

  if (sequence.length > maxSuitableLength) {
    hasProblem = true;
    report = report + "Contains more than " + maxSuitableLength + " bases; ";
  }

  if (hasProblem) {
    return "Warning: " + report;
  } else {
    return "Pass";
  }
}

function _getSuitableTmReport(sequence, range, minSuitable, maxSuitable) {
  var report = "";
  var hasProblem = false;
  var lowerCalculated;
  var upperCalculated;

  //range may be a single number or a string containing something like "40 to 60";
  var rangePattern = new RegExp("([d.]+)D+([d.]+)", "gi");
  if (range.search(rangePattern) != -1) {
    lowerCalculated = parseFloat($1);
    upperCalculated = parseFloat($2);
  } else {
    lowerCalculated = parseFloat(range);
    upperCalculated = parseFloat(range);
  }

  if (lowerCalculated < minSuitable) {
    hasProblem = true;
    report = report + "Tm is less than " + minSuitable + "; ";
  }

  if (upperCalculated > maxSuitable) {
    hasProblem = true;
    report = report + "Tm is greater than " + maxSuitable + "; ";
  }

  if (hasProblem) {
    return "Warning: " + report;
  } else {
    return "Pass";
  }
}

function _getSuitableThreePrimeGC(sequence, minSuitable, maxSuitable) {
  var threePrimeEnd;
  var desiredThreePrimeEndSize = 5;
  var report = "";
  var gcCounts;
  var hasProblem = false;
  if (sequence.length >= desiredThreePrimeEndSize) {
    threePrimeEnd = sequence.substr(
      sequence.length - desiredThreePrimeEndSize,
      5
    );
  } else {
    threePrimeEnd = sequence;
  }

  gcCounts =
    _getBaseCount(threePrimeEnd, "g") + _getBaseCount(threePrimeEnd, "c");

  if (gcCounts < minSuitable) {
    hasProblem = true;
    report =
      report +
      "There are less than " +
      minSuitable +
      " G's or C's in the last " +
      threePrimeEnd.length +
      " bases; ";
  }

  if (gcCounts > maxSuitable) {
    hasProblem = true;
    report =
      report +
      "There are more than " +
      maxSuitable +
      " G's or C's in the last " +
      threePrimeEnd.length +
      " bases; ";
  }

  if (hasProblem) {
    return "Warning: " + report;
  } else {
    return "Pass";
  }
}

function _getSelfComplementarityReport(sequence, maxContig, maxPercentIdent) {
  var matchScore = 1;
  var mismatchScore = -1;
  var gapPenalty = 99;
  var beginGapPenalty = 0;
  var endGapPenalty = 0;

  var returnHash = {};
  returnHash["report"] = "";
  returnHash["upper"] = "";
  returnHash["lower"] = "";
  returnHash["divider"] = "";

  var report = "";
  var hasProblem = false;
  var sequenceLength = sequence.length;

  var matrix = new Complementarity();
  matrix.setMatch(matchScore);
  matrix.setMismatch(mismatchScore);

  var scoreSet = new ScoreSet();
  scoreSet.setScoreSetParam(matrix, gapPenalty, beginGapPenalty, endGapPenalty);

  var rev = reverse(sequence);
  //convert String to Array
  sequence = sequence.match(/./g);
  rev = rev.match(/./g);

  //align_pair_quad.js
  alignment = new AlignPairQuad();
  alignment.initializeMatrix(sequence, rev, scoreSet);
  alignment.fillMatrix();
  alignment.align();

  //aligned output will be something like:
  //cttttgagcaagttcagcctggttaag--
  //--gaattggtccgacttgaacgagttttc
  var seqAligned = alignment.getAlignedM().replace(/\-/g, " ");
  var revAligned = alignment.getAlignedN().replace(/\-/g, " ");

  var score = alignment.score;

  var divider = new Array();
  var maxContiguous = 0;
  var totalMatches = 0;
  var contiguous = 0;
  for (var i = 0; i < seqAligned.length; i++) {
    if (
      scoreSet.getScore(seqAligned.charAt(i), revAligned.charAt(i)) ==
      matchScore
    ) {
      divider.push("|");
      contiguous++;
      totalMatches++;
    } else {
      divider.push(" ");
      contiguous = 0;
    }

    if (contiguous > maxContiguous) {
      maxContiguous = contiguous;
    }
  }

  if (maxContiguous > maxContig) {
    hasProblem = true;
    report =
      report +
      "There are more than " +
      maxContig +
      " self-annealing bases in a row; ";
  }

  if ((totalMatches / sequenceLength) * 100 > maxPercentIdent) {
    hasProblem = true;
    report =
      report +
      "More than " +
      maxPercentIdent +
      " percent of the bases are self-annealing; ";
  }

  if (hasProblem) {
    report = "Warning: " + report;
  } else {
    report = "Pass";
  }

  returnHash["report"] = report;
  returnHash["upper"] = seqAligned;
  returnHash["lower"] = revAligned;
  returnHash["divider"] = divider.join("");

  return returnHash;
}

function _getHairpinReport(sequence, maxContig, maxPercentIdent) {
  var upper = sequence;
  upper = upper.match(/./g);
  var loop = "";
  var lower = new Array();

  var returnHash = {};
  returnHash["report"] = "";
  returnHash["upper"] = "";
  returnHash["lower"] = "";
  returnHash["divider"] = "";

  var topScore = 0;
  var score;
  var u;
  var l;
  var topScoreUpper = sequence;
  var topScoreLower = "";
  var topScoreLoop = "";

  var matchScore = 1;
  var mismatchScore = -1;
  var gapPenalty = 99;
  var beginGapPenalty = 0;
  var endGapPenalty = 0;

  var report = "";
  var hasProblem = false;
  var sequenceLength = sequence.length;

  var matrix = new Complementarity();
  matrix.setMatch(matchScore);
  matrix.setMismatch(mismatchScore);

  var scoreSet = new ScoreSet();
  scoreSet.setScoreSetParam(matrix, gapPenalty, beginGapPenalty, endGapPenalty);

  while (upper.length > 0) {
    score = 0;
    if (loop == "") {
      loop = upper.pop();
    } else {
      lower.push(loop);
      loop = "";
    }

    //determine score
    u = upper.length - 1;
    l = lower.length - 1;
    while (u >= 0 && l >= 0) {
      score = score + scoreSet.getScore(upper[u], lower[l]);
      u--;
      l--;
    }

    if (score > topScore) {
      topScore = score;
      topScoreUpper = upper.join("");
      topScoreLower = lower.join("");
      topScoreLoop = loop;
    }
  }

  //format top scoring hit and return

  var upperLowerDiff = topScoreUpper.length - topScoreLower.length;
  if (upperLowerDiff > 0) {
    for (var i = 0; i < upperLowerDiff; i++) {
      topScoreLower = " " + topScoreLower;
    }
  } else if (upperLowerDiff < 0) {
    for (var i = upperLowerDiff; i < 0; i++) {
      topScoreUpper = " " + topScoreUpper;
    }
  }

  var maxContiguous = 0;
  var totalMatches = 0;
  var contiguous = 0;
  var divider = new Array();
  //add vertical lines between matches
  for (var i = 0; i < topScoreUpper.length; i++) {
    if (
      scoreSet.getScore(topScoreUpper.charAt(i), topScoreLower.charAt(i)) ==
      matchScore
    ) {
      divider.push("|");
      contiguous++;
      totalMatches++;
    } else {
      divider.push(" ");
      contiguous = 0;
    }

    if (contiguous > maxContiguous) {
      maxContiguous = contiguous;
    }
  }

  if (maxContiguous > maxContig) {
    hasProblem = true;
    report =
      report + "There are more than " + maxContig + " hairpin bases in a row; ";
  }

  if ((totalMatches / sequenceLength) * 100 > maxPercentIdent) {
    hasProblem = true;
    report =
      report +
      "More than " +
      maxPercentIdent +
      " percent of the bases are in a hairpin; ";
  }

  if (hasProblem) {
    report = "Warning: " + report;
  } else {
    report = "Pass";
  }

  if (topScoreLoop == "") {
    topScoreLoop = ")";
  }

  returnHash["report"] = report;
  returnHash["upper"] = topScoreUpper;
  returnHash["divider"] = divider.join("") + topScoreLoop;
  returnHash["lower"] = topScoreLower;

  return returnHash;
}

//I wrote these classes originally in pairwise_align_dna.js
//------------------------------------ ScoreSet class

//ScoreSet getScore
function getScore(r1, r2) {
  return this.scoringMatrix.scoringMatrix_getScore(r1, r2);
}

//ScoreSet setScoreSetParam
function setScoreSetParam(
  scoringMatrix,
  gapPenalty,
  beginGapPenalty,
  endGapPenalty
) {
  this.scoringMatrix = scoringMatrix;
  this.gap = gapPenalty;
  this.beginGap = beginGapPenalty;
  this.endGap = endGapPenalty;
}

//ScoreSet class
function ScoreSet() {
  this.scoringMatrix;
  this.gap;
  this.beginGap;
  this.endGap;
  this.useBeginGapTop = true;
  this.useBeginGapLeft = true;
  this.useEndGapBottom = true;
  this.useEndGapRight = true;
}

//create and throw away a prototype object
new ScoreSet();

//define object methods
ScoreSet.prototype.getScore = getScore;
ScoreSet.prototype.setScoreSetParam = setScoreSetParam;

//------------------------------------

//------------------------------------ ScoringMatrix abstract class
//ScoringMatrix getScore method
function scoringMatrix_getScore(r1, r2) {
  r1 = r1.toLowerCase();
  r2 = r2.toLowerCase();
  if ((r1 == "g" && r2 == "c") || (r2 == "g" && r1 == "c")) {
    return this.match;
  } else if ((r1 == "a" && r2 == "t") || (r2 == "a" && r1 == "t")) {
    return this.match;
  } else {
    return this.mismatch;
  }
}

//ScoringMatrix class
function ScoringMatrix() {
  this.mismatch;
  this.match;
}

//create and throw away a prototype object
new ScoringMatrix();

//define object methods
ScoringMatrix.prototype.scoringMatrix_getScore = scoringMatrix_getScore;

//------------------------------------ Complementarity class extends ScoringMatrix Class
//Complementarity class setMismatch method
function setMismatch(mismatchScore) {
  this.mismatch = mismatchScore;
}

//Complementarity class setMatch method
function setMatch(matchScore) {
  this.match = matchScore;
}

//Complementarity class
function Complementarity() {}

Complementarity.prototype = new ScoringMatrix();
Complementarity.prototype.setMismatch = setMismatch;
Complementarity.prototype.setMatch = setMatch;
