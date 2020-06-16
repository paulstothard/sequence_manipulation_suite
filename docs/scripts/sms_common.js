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

function addReturns(sequence) {
  sequence = sequence.replace(/(.{60})/g, function (str, p1, offset, s) {
    return p1 + "\n";
  });
  return sequence;
}

function checkAlign(arrayOfTitles, arrayOfSequences) {
  var lengthOfAlign = arrayOfSequences[0].length;
  if (arrayOfSequences.length < 2) {
    alert("Please enter an alignment consisting of at least two sequences.");
    return false;
  }
  for (var i = 0; i < arrayOfTitles.length; i++) {
    if (
      arrayOfTitles[i].search(/\S/) == -1 ||
      arrayOfSequences[i].search(/\S/) == -1 ||
      arrayOfSequences[i].length != lengthOfAlign
    ) {
      alert("There is a problem with the alignment format.");
      return false;
    }
  }
  return true;
}

function checkCodonTable(codonTable) {
  if (
    codonTable.search(/AmAcid/) == -1 ||
    codonTable.search(/Codon/) == -1 ||
    codonTable.search(/Number/) == -1 ||
    codonTable.search(/\/1000/) == -1 ||
    codonTable.search(/Fraction\s*\.\./) == -1
  ) {
    alert("The codon table has been entered incorrectly.");
    return false;
  }
  return true;
}

function checkFormElement(formElement) {
  if (formElement.value.search(/\S/) == -1) {
    alert("Please enter some text.");
    return false;
  }
  return true;
}

function checkGeneticCode(arrayOfPatterns) {
  var z = 0;
  var codon = "";
  var oneMatch = false;
  var testSequence =
    "gggggaggtggcgaggaagatgacgtggtagttgtcgcggcagctgccaggagaagtagcaagaaaaataacatgataattatcacgacaactacctggtgatgttgctagtaatattacttgttatttttctcgtcatcttcccggcgacgtcgccagcaacatcacctgctacttctcccgccacctccc";
  while (z < arrayOfPatterns.length) {
    if (arrayOfPatterns[z].search(/^\s*\/[a-zA-Z\|\[\]]+\/=[a-zA-Z\*]/) == -1) {
      alert(
        "Genetic code error: one or more patterns have been entered incorrectly."
      );
      return false;
    }
    if (moreExpressionCheck(arrayOfPatterns[z]) == false) {
      alert(
        "Genetic code error: one or more patterns have been entered incorrectly."
      );
      return false;
    }
    z = z + 1;
  }
  var geneticCodeMatchResult = new Array(arrayOfPatterns.length);
  var geneticCodeMatchExp = new Array(arrayOfPatterns.length);
  for (var j = 0; j < arrayOfPatterns.length; j++) {
    geneticCodeMatchExp[j] = eval(arrayOfPatterns[j].match(/\/.+\//) + "gi");
    geneticCodeMatchResult[j] = arrayOfPatterns[j]
      .match(/=[a-zA-Z\*]/)
      .toString();
    geneticCodeMatchResult[j] = geneticCodeMatchResult[j].replace(/=/g, "");
  }
  for (var i = 0; i <= testSequence.length - 3; i = i + 3) {
    codon = testSequence.substring(i, i + 3);
    for (var j = 0; j < geneticCodeMatchExp.length; j++) {
      if (codon.search(geneticCodeMatchExp[j]) != -1) {
        if (oneMatch == true) {
          alert(
            "Genetic code error: more than one amino acid is coded by the codon: " +
              codon +
              "."
          );
          return false;
        }
        oneMatch = true;
      }
    }
    if (oneMatch == false) {
      alert("The genetic code expressions are missing a codon.");
      return false;
    }
    oneMatch = false;
  }
  return true;
}

function checkGroupInput(arrayOfPatterns) {
  var z = 0;
  while (z < arrayOfPatterns.length) {
    if (arrayOfPatterns[z].search(/[^acdefghiklmnpqrstvwyz]/i) != -1) {
      alert("One or more groups have been entered incorrectly.");
      return false;
    }
    z = z + 1;
  }
  for (var i = 0; i < arrayOfPatterns.length; i++) {
    var re = new RegExp("[" + arrayOfPatterns[i] + "]", "gi");
    for (var j = i + 1; j < arrayOfPatterns.length; j++) {
      if (arrayOfPatterns[j].search(re) != -1) {
        alert("The same amino acid is in more than one similarity group.");
        return false;
      }
    }
  }
  return true;
}

function checkRestPatterns(arrayOfPatterns) {
  var z = 0;
  while (z < arrayOfPatterns.length) {
    if (
      arrayOfPatterns[z].search(
        /^\s*\/[acgturyswkmbdhvn\[\]]+\/\s+\([^\/]+\)\d+/i
      ) == -1
    ) {
      alert("One or more patterns have been entered incorrectly.");
      return false;
    }
    if (moreExpressionCheck(arrayOfPatterns[z]) == false) {
      alert("One or more patterns have been entered incorrectly.");
      return false;
    }
    z = z + 1;
  }
  return true;
}

function checkSequenceLength(text, maxInput) {
  if (getSequenceFromFasta(text).replace(/[^A-Za-z]/g, "").length > maxInput) {
    alert(
      "Please enter a sequence consisting of less than or equal to " +
        maxInput +
        " characters."
    );
    return false;
  } else {
    return true;
  }
}

function checkTextLength(text, maxInput) {
  if (text.length > maxInput) {
    alert(
      "Please enter text consisting of less than or equal to " +
        maxInput +
        " characters."
    );
    return false;
  } else {
    return true;
  }
}

function complement(dnaSequence) {
  //there is no tr operator
  //should write a tr method to replace this
  dnaSequence = dnaSequence.replace(/g/g, "1");
  dnaSequence = dnaSequence.replace(/c/g, "2");
  dnaSequence = dnaSequence.replace(/1/g, "c");
  dnaSequence = dnaSequence.replace(/2/g, "g");
  dnaSequence = dnaSequence.replace(/G/g, "1");
  dnaSequence = dnaSequence.replace(/C/g, "2");
  dnaSequence = dnaSequence.replace(/1/g, "C");
  dnaSequence = dnaSequence.replace(/2/g, "G");

  dnaSequence = dnaSequence.replace(/a/g, "1");
  dnaSequence = dnaSequence.replace(/t/g, "2");
  dnaSequence = dnaSequence.replace(/1/g, "t");
  dnaSequence = dnaSequence.replace(/2/g, "a");
  dnaSequence = dnaSequence.replace(/A/g, "1");
  dnaSequence = dnaSequence.replace(/T/g, "2");
  dnaSequence = dnaSequence.replace(/1/g, "T");
  dnaSequence = dnaSequence.replace(/2/g, "A");

  dnaSequence = dnaSequence.replace(/u/g, "a");
  dnaSequence = dnaSequence.replace(/U/g, "A");

  dnaSequence = dnaSequence.replace(/r/g, "1");
  dnaSequence = dnaSequence.replace(/y/g, "2");
  dnaSequence = dnaSequence.replace(/1/g, "y");
  dnaSequence = dnaSequence.replace(/2/g, "r");
  dnaSequence = dnaSequence.replace(/R/g, "1");
  dnaSequence = dnaSequence.replace(/Y/g, "2");
  dnaSequence = dnaSequence.replace(/1/g, "Y");
  dnaSequence = dnaSequence.replace(/2/g, "R");

  dnaSequence = dnaSequence.replace(/k/g, "1");
  dnaSequence = dnaSequence.replace(/m/g, "2");
  dnaSequence = dnaSequence.replace(/1/g, "m");
  dnaSequence = dnaSequence.replace(/2/g, "k");
  dnaSequence = dnaSequence.replace(/K/g, "1");
  dnaSequence = dnaSequence.replace(/M/g, "2");
  dnaSequence = dnaSequence.replace(/1/g, "M");
  dnaSequence = dnaSequence.replace(/2/g, "K");

  dnaSequence = dnaSequence.replace(/b/g, "1");
  dnaSequence = dnaSequence.replace(/v/g, "2");
  dnaSequence = dnaSequence.replace(/1/g, "v");
  dnaSequence = dnaSequence.replace(/2/g, "b");
  dnaSequence = dnaSequence.replace(/B/g, "1");
  dnaSequence = dnaSequence.replace(/V/g, "2");
  dnaSequence = dnaSequence.replace(/1/g, "V");
  dnaSequence = dnaSequence.replace(/2/g, "B");

  dnaSequence = dnaSequence.replace(/d/g, "1");
  dnaSequence = dnaSequence.replace(/h/g, "2");
  dnaSequence = dnaSequence.replace(/1/g, "h");
  dnaSequence = dnaSequence.replace(/2/g, "d");
  dnaSequence = dnaSequence.replace(/D/g, "1");
  dnaSequence = dnaSequence.replace(/H/g, "2");
  dnaSequence = dnaSequence.replace(/1/g, "H");
  dnaSequence = dnaSequence.replace(/2/g, "D");

  return dnaSequence;
}

function closeForm() {
  outputWindow.document.write("</form>");
  return true;
}

function closePre() {
  //this is needed for Internet Explorer bug
  outputWindow.document.write("</div>");
  outputWindow.document.write("</pre>\n");
}

function closeTextArea() {
  outputWindow.document.write("</textarea>");
  return true;
}

function closeWindow() {
  outputWindow.document.write("</body>\n</html>\n");
  outputWindow.status = "Done.";
  outputWindow.document.close();
  return true;
}

function convertDegenerates(sequence) {
  sequence = sequence.toLowerCase();
  sequence = sequence.replace(/t/g, "[TU]");
  sequence = sequence.replace(/r/g, "[AGR]");
  sequence = sequence.replace(/y/g, "[CTUY]");
  sequence = sequence.replace(/s/g, "[GCS]");
  sequence = sequence.replace(/w/g, "[ATUW]");
  sequence = sequence.replace(/k/g, "[GTUK]");
  sequence = sequence.replace(/m/g, "[ACM]");
  sequence = sequence.replace(/b/g, "[CGTUBSKY]");
  sequence = sequence.replace(/d/g, "[AGTUDRKW]");
  sequence = sequence.replace(/h/g, "[ACTUHMYW]");
  sequence = sequence.replace(/v/g, "[ACGVSMR]");
  sequence = sequence.replace(/n/g, "[ACGTURYSWKMBDHVN]");
  return sequence;
}

function earlyCheckAlign(alignArray) {
  if (alignArray.length < 3) {
    alert("There is a problem with the alignment format.");
    return false;
  }
  for (var i = 1; i < alignArray.length; i++) {
    if (alignArray[i].search(/[^\s]+\s/) == -1) {
      alert("There is a problem with the alignment format.");
      return false;
    }
  }
  return true;
}

function filterAlignSeq(alignSeq) {
  alignSeq = alignSeq.replace(/[^abcdefghiklmnpqrstvwxyz\.\-]/gi, "");
  return alignSeq;
}

function filterAlignSeqAllowAsterisk(alignSeq) {
  alignSeq = alignSeq.replace(/[^abcdefghiklmnpqrstvwxyz\.\-\*]/gi, "");
  return alignSeq;
}

function filterFastaTitle(sequenceTitle) {
  sequenceTitle = sequenceTitle.replace(/\s{2,}/g, " ");
  sequenceTitle = sequenceTitle.replace(/^\s*/g, "");
  sequenceTitle = sequenceTitle.replace(/[\f\n\r\t]+$/g, "\n");
  return sequenceTitle.replace(/[\<\>]\n/gi, "");
}

function getArrayOfFasta(sequenceData) {
  var arrayOfFasta = new Array();
  var matchArray;
  var re = /\>[^\>]+/g;
  if (sequenceData.search(/\>[^\f\n\r]+[\f\n\r]/) != -1) {
    while ((matchArray = re.exec(sequenceData))) {
      arrayOfFasta.push(matchArray[0]);
    }
  } else {
    arrayOfFasta[0] = sequenceData;
  }
  return arrayOfFasta;
}

function getFastaTitleFromTitleAndSequence(fastaSequenceTitle, sequence) {
  var stringToReturn =
    "&gt;results for " + sequence.length + " residue sequence ";
  if (fastaSequenceTitle.search(/[^\s]/) != -1) {
    stringToReturn = stringToReturn + '"' + fastaSequenceTitle + '"';
  }
  stringToReturn =
    stringToReturn + ' starting "' + sequence.substring(0, 10) + '"';
  return stringToReturn + "\n";
}

function getFuzzySearchTitle(
  fastaSequenceTitleOne,
  sequenceOne,
  fastaSequenceTitleTwo,
  sequenceTwo
) {
  var stringToReturn =
    "Search results for " + sequenceOne.length + " residue sequence ";
  if (fastaSequenceTitleOne.search(/[^\s]/) != -1) {
    stringToReturn = stringToReturn + '"' + fastaSequenceTitleOne + '"';
  }
  stringToReturn =
    stringToReturn + ' starting "' + sequenceOne.substring(0, 10) + '"\n';
  stringToReturn =
    stringToReturn + "and " + sequenceTwo.length + " residue sequence ";
  if (fastaSequenceTitleTwo.search(/[^\s]/) != -1) {
    stringToReturn = stringToReturn + '"' + fastaSequenceTitleTwo + '"';
  }
  stringToReturn =
    stringToReturn + ' starting "' + sequenceTwo.substring(0, 10) + '"';
  return '<div class="info">' + stringToReturn + "</div>\n";
}

function getGeneticCodeMatchExp(arrayOfPatterns) {
  var geneticCodeMatchExp = new Array(arrayOfPatterns.length);
  for (var j = 0; j < arrayOfPatterns.length; j++) {
    geneticCodeMatchExp[j] = eval(arrayOfPatterns[j].match(/\/.+\//) + "gi");
  }
  return geneticCodeMatchExp;
}

function getGeneticCodeMatchResult(arrayOfPatterns) {
  var geneticCodeMatchResult = new Array(arrayOfPatterns.length);
  for (var j = 0; j < arrayOfPatterns.length; j++) {
    geneticCodeMatchResult[j] = arrayOfPatterns[j]
      .match(/=[a-zA-Z\*]/)
      .toString();
    geneticCodeMatchResult[j] = geneticCodeMatchResult[j].replace(/=/g, "");
  }
  return geneticCodeMatchResult;
}

function getInfoFromTitleAndSequence(fastaSequenceTitle, sequence) {
  var stringToReturn = "Results for " + sequence.length + " residue sequence ";
  if (fastaSequenceTitle.search(/[^\s]/) != -1) {
    stringToReturn = stringToReturn + '"' + fastaSequenceTitle + '"';
  }
  stringToReturn =
    stringToReturn + ' starting "' + sequence.substring(0, 10) + '"';
  return '<div class="info">' + stringToReturn + "</div>\n";
}

function getInfoFromTitleAndSequenceAndTopology(
  fastaSequenceTitle,
  sequence,
  topology
) {
  var stringToReturn =
    "Results for " + topology + " " + sequence.length + " residue sequence ";
  if (fastaSequenceTitle.search(/[^\s]/) != -1) {
    stringToReturn = stringToReturn + '"' + fastaSequenceTitle + '"';
  }
  stringToReturn =
    stringToReturn + ' starting "' + sequence.substring(0, 10) + '"';
  return '<div class="info">' + stringToReturn + "</div>\n";
}

function getPairwiseAlignTitle(
  fastaSequenceTitleOne,
  sequenceOne,
  fastaSequenceTitleTwo,
  sequenceTwo
) {
  var stringToReturn =
    "Alignment results for " + sequenceOne.length + " residue sequence ";
  if (fastaSequenceTitleOne.search(/[^\s]/) != -1) {
    stringToReturn = stringToReturn + '"' + fastaSequenceTitleOne + '"';
  }
  stringToReturn =
    stringToReturn + ' starting "' + sequenceOne.substring(0, 10) + '"\n';
  stringToReturn =
    stringToReturn + "and " + sequenceTwo.length + " residue sequence ";
  if (fastaSequenceTitleTwo.search(/[^\s]/) != -1) {
    stringToReturn = stringToReturn + '"' + fastaSequenceTitleTwo + '"';
  }
  stringToReturn =
    stringToReturn + ' starting "' + sequenceTwo.substring(0, 10) + '"';
  return '<div class="info">' + stringToReturn + "</div>\n";
}

function getRandomSequence(components, lengthOut) {
  var sequenceArray = new Array();
  var tempNum = 0;
  var tempChar = "";
  for (var j = 0; j < lengthOut; j++) {
    tempNum = Math.floor(Math.random() * components.length);
    tempChar = components[tempNum];
    sequenceArray.push(tempChar);
  }
  return sequenceArray.join("");
}

function getSequenceFromFasta(sequenceRecord) {
  if (sequenceRecord.search(/\>[^\f\n\r]+[\f\n\r]/) != -1) {
    sequenceRecord = sequenceRecord.replace(/\>[^\f\n\r]+[\f\n\r]/, "");
  }
  return sequenceRecord;
}

function getTitleFromFasta(sequenceRecord) {
  var fastaTitle = "Untitled";
  if (sequenceRecord.search(/\>[^\f\n\r]+[\f\n\r]/) != -1) {
    fastaTitle = sequenceRecord.match(/\>[^\f\n\r]+[\f\n\r]/, "").toString();
    fastaTitle = fastaTitle.replace(/\>|[\f\n\r]/g, "");
    fastaTitle = fastaTitle.replace(/\s{2,}/g, " ");
    fastaTitle = fastaTitle.replace(/[\<\>]/gi, "");
  }
  return fastaTitle;
}

function moreExpressionCheck(expressionToCheck) {
  if (
    expressionToCheck.search(/\[[A-Za-z\|]*\[/) != -1 ||
    expressionToCheck.search(/\][A-Za-z\|]*\]/) != -1 ||
    expressionToCheck.search(/\[\]/) != -1 ||
    expressionToCheck.search(/\/[A-Za-z\|]*\]/) != -1 ||
    expressionToCheck.search(/\[[A-Za-z\|]*\//) != -1 ||
    expressionToCheck.search(/\|\|/) != -1 ||
    expressionToCheck.search(/\/\|/) != -1 ||
    expressionToCheck.search(/\|\//) != -1 ||
    expressionToCheck.search(/\[.\]/) != -1 ||
    expressionToCheck.search(/\</) != -1 ||
    expressionToCheck.search(/\>/) != -1
  ) {
    return false;
  }
  return true;
}

function openForm() {
  outputWindow.document.write('<form action="">\n');
  return true;
}

function openPre() {
  //this is needed for Internet Explorer bug
  outputWindow.document.write("<pre>");
  outputWindow.document.write('<div class="pre">');
}

function openTextArea() {
  outputWindow.document.write('<br /><textarea rows="6" cols="61">\n');
  return true;
}

function openWindow(title) {
  _openWindow(title, true);
}

function _openWindow(title, isColor) {
  outputWindow = window.open(
    "",
    "my_new_window",
    "toolbar=no, location=no, directories=no, status=yes, menubar=yes, scrollbars=yes, resizable=yes, copyhistory=no, width=800, height=400"
  );
  outputWindow.focus();
  outputWindow.document.write(
    '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">\n' +
      '<html lang="en">\n' +
      "<head>\n" +
      "<title>Sequence Manipulation Suite</title>\n" +
      '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1" />\n'
  );

  if (isColor) {
    outputWindow.document.write(
      '<style type="text/css">\n' +
        "body.main {font-size:90%; font-family: arial, sans-serif; color: #000000; background-color: #FFFFFF}\n" +
        "div.pre {color: #000000; font-family: courier, sans-serif; white-space: pre}\n" +
        "div.title {color: #000000; text-align: left; background-color: #FFFFFF}\n" +
        "div.info {font-weight: bold}\n" +
        "span.none, td.none {color: #000000; background-color: #FFFFFF}\n" +
        "span.one, td.one {color: #000000; background-color: #66FF00}\n" +
        "span.two, td.two {color: #000000; background-color: #FFFF66}\n" +
        "span.three, td.three {color: #000000; background-color: #FFFFFF}\n" +
        "span.forward_primer, td.forward_primer {color: #000000; background-color: #FF66FF}\n" +
        "span.reverse_primer, td.reverse_primer {color: #000000; background-color: #FF9933}\n" +
        "span.current_sequence {color: #000000; background-color: #FFFFFF}\n" +
        "span.mutated_sequence {color: #990066; background-color: #FFFFFF}\n" +
        "td.many {color: #000000}\n" +
        "td.title {font-weight: bold; color: #000000; background-color: #FFFFFF}\n" +
        "</style>\n"
    );
  } else {
    outputWindow.document.write(
      '<style type="text/css">\n' +
        "body.main {font-size:90%; font-family: arial, sans-serif; color: #000000; background-color: #FFFFFF; margin: 0 auto; padding: 0}\n" +
        "div.pre {color: #000000; background-color: #FFFFFF; font-family: courier, sans-serif; white-space: pre}\n" +
        "div.title {display: none}\n" +
        "div.info {font-weight: bold}\n" +
        "span.none, td.none {color: #000000; background-color: #FFFFFF}\n" +
        "span.one, td.one {color: #000000; text-decoration: underline; background-color: #FFFFFF}\n" +
        "span.two, td.two {color: #000000; font-style: italic; background-color: #FFFFFF}\n" +
        "span.three, td.three {color: #000000; background-color: #FFFFFF}\n" +
        "span.forward_primer, td.forward_primer {color: #000000; background-color: #FFFFFF}\n" +
        "span.reverse_primer, td.reverse_primer {color: #000000; background-color: #FFFFFF}\n" +
        "span.current_sequence {color: #000000; background-color: #FFFFFF}\n" +
        "span.mutated_sequence {color: #000000; text-decoration: underline; background-color: #FFFFFF}\n" +
        "td.many {color: #000000; background-color: #FFFFFF}\n" +
        "td.title {font-weight: bold; color: #000000; background-color: #FFFFFF}\n" +
        "img {display: none}\n" +
        "</style>\n"
    );
  }
  outputWindow.document.write(
    "</head>\n" +
      '<body class="main">\n' +
      '<div class="title">' +
      title +
      " results</div>\n"
  );
  outputWindow.status = "Please Wait.";
  return true;
}

function openWindowAlign(title) {
  _openWindowAlign(title, true);
}

function _openWindowAlign(title, isBackground) {
  outputWindow = window.open(
    "",
    "my_new_window",
    "toolbar=no, location=no, directories=no, status=yes, menubar=yes, scrollbars=yes, resizable=yes, copyhistory=no, width=800, height=400"
  );
  outputWindow.focus();
  outputWindow.document.write(
    '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">\n' +
      '<html lang="en">\n' +
      "<head>\n" +
      "<title>Sequence Manipulation Suite</title>\n" +
      '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1" />\n'
  );
  if (isBackground) {
    outputWindow.document.write(
      '<style type="text/css">\n' +
        "body.main {font-family: arial, sans-serif; font-size:90%; color: #000000; background-color: #FFFFFF}\n" +
        "div.pre {color: #000000; font-family: courier, sans-serif; white-space: pre}\n" +
        "div.title {color: #000000; text-align: left; background-color: #FFFFFF}\n" +
        "div.info {font-weight: bold}\n" +
        "span.ident {color: #FFFFFF; background-color: #000000}\n" +
        "span.sim {color: #FFFFFF; background-color: #666666}\n" +
        "span.g, span.a, span.v, span.l, span.i {color: #000000; background-color: #C0C0C0}\n" +
        "span.f, span.y, span.w {color: #000000; background-color: #FF6600}\n" +
        "span.c, span.m {color: #000000; background-color: #FFFF00}\n" +
        "span.s, span.t {color: #000000; background-color: #66FF00}\n" +
        "span.k, span.r, span.h {color: #000000; background-color: #FF0000}\n" +
        "span.d, span.e {color: #000000; background-color: #0066FF}\n" +
        "span.n, span.q {color: #000000; background-color: #996633}\n" +
        "span.p {color: #000000; background-color: #FF99FF}\n" +
        "</style>\n"
    );
  } else {
    outputWindow.document.write(
      '<style type="text/css">\n' +
        "body.main {font-family: arial, sans-serif; font-size:90%; color: #000000; background-color: #FFFFFF}\n" +
        "div.pre {color: #000000; font-family: courier, sans-serif; white-space: pre}\n" +
        "div.title {display: none}\n" +
        "div.info {font-weight: bold}\n" +
        "span.ident {color: #000000; font-weight: bold; text-decoration: underline; background-color: #FFFFFF}\n" +
        "span.sim {color: #000000; font-weight: bold; background-color: #FFFFFF}\n" +
        "span.diff {color: #999999; background-color: #FFFFFF}\n" +
        "span.g, span.a, span.v, span.l, span.i {color: #CC33CC; background-color: #FFFFFF}\n" +
        "span.f, span.y, span.w {color: #FF6600; background-color: #FFFFFF}\n" +
        "span.c, span.m {color: #FFCC00; background-color: #FFFFFF}\n" +
        "span.s, span.t {color: #CCFF00; background-color: #FFFFFF}\n" +
        "span.k, span.r, span.h {color: #FF0000; background-color: #FFFFFF}\n" +
        "span.d, span.e {color: #0000FF; background-color: #FFFFFF}\n" +
        "span.n, span.q {color: #996633; background-color: #FFFFFF}\n" +
        "span.p {color: #00FFCC; background-color: #FFFFFF}\n" +
        "img {display: none}\n" +
        "</style>\n"
    );
  }
  outputWindow.document.write(
    "</head>\n" +
      '<body class="main">\n' +
      '<div class="title">' +
      title +
      " results</div>\n"
  );
  outputWindow.status = "Please Wait.";
  return true;
}

function removeFormatting(sequence) {
  return sequence.replace(/[\d\s]/g, "");
}

function removeNonDna(sequence) {
  return sequence.replace(/[^gatucryswkmbdhvnxGATUCRYSWKMBDHVNX]/g, "");
}

function removeNonDnaStrict(sequence) {
  return sequence.replace(/[^gatucGATUC]/g, "");
}

function removeNonProtein(sequence) {
  return sequence.replace(
    /[^ACDEFGHIKLMNPQRSTVWYZacdefghiklmnpqrstvwyz\*]/g,
    ""
  );
}

function removeNonProteinStrict(sequence) {
  return sequence.replace(
    /[^ACDEFGHIKLMNPQRSTVWYZacdefghiklmnpqrstvwyz\*]/g,
    ""
  );
}

function removeNonProteinAllowDegen(sequence) {
  return sequence.replace(
    /[^ABCDEFGHIKLMNPQRSTVWYXZabcdefghiklmnpqrstvwyxz\*]/g,
    ""
  );
}

function removeNonProteinAllowX(sequence) {
  return sequence.replace(
    /[^ACDEFGHIKLMNPQRSTVWYZXacdefghiklmnpqrstvwyzx\*]/g,
    ""
  );
}

function removeWhiteSpace(text) {
  return text.replace(/\s/g, "");
}

function removeNonLetters(sequence) {
  return sequence.replace(/[^A-Z]/gi, "");
}

function reverse(dnaSequence) {
  var tempDnaArray = new Array();
  if (dnaSequence.search(/./) != -1) {
    tempDnaArray = dnaSequence.match(/./g);
    tempDnaArray = tempDnaArray.reverse();
    dnaSequence = tempDnaArray.join("");
  }
  return dnaSequence;
}

function rightNum(theNumber, sequenceToAppend, lengthOfColumn, tabIn) {
  var j = 0;
  var tempString = "";
  theNumber = theNumber.toString();
  for (var j = theNumber.length; j < lengthOfColumn; j++) {
    tempString = tempString + " ";
  }
  theNumber = tempString + theNumber + " ";
  sequenceToAppend = sequenceToAppend + theNumber + tabIn;
  return sequenceToAppend;
}

function testScript() {
  //test some javascript functions to see how the browser performs.
  //want to prevent non Javascript 1.5 browsers from attempting to run.
  //first test Array.push()
  var testArray = new Array();
  var testString = "1234567890";

  testArray.push(testString);
  if (testArray[0] != testString) {
    alert(
      "Array object push method not supported. See browser compatibility page."
    );
    return false;
  }

  //now test the 'm' flag in a regular expression
  testString = "1\n2\n3";
  var re = /^2$/m;
  if (testString.search(re) == -1) {
    alert(
      "Regular expression 'm' flag not supported. See browser compatibility page."
    );
    return false;
  }

  var caughtException = false;
  //now test exception handling
  try {
    re = eval(
      "Exception handling not supported. Check browser compatibility page."
    );
  } catch (e) {
    caughtException = true;
  }

  if (!caughtException) {
    alert("Exception handling not supported. See browser compatibility page.");
  }

  //now test replace lambda function
  testString = "123";
  testString = testString.replace(/(\d)/g, function (str, p1, offset, s) {
    return p1 + "X";
  });
  if (testString != "1X2X3X") {
    alert(
      "Nested function in String replace method not supported. See browser compatibility page."
    );
    return false;
  }

  //test number methods toFixed() and toPrecision()
  var testNum = 2489.8237;
  if (testNum.toFixed(3) != 2489.824) {
    alert(
      "Number toFixed() method not supported. See browser compatibility page."
    );
    return false;
  }

  if (testNum.toPrecision(5) != 2489.8) {
    alert(
      "Number toPrecision() method not supported. See browser compatibility page."
    );
    return false;
  }

  return true;
}

function verifyDigits(theNumber) {
  if (theNumber.search(/\d/) == -1) {
    alert("Please enter a number");
    return false;
  }
}

function verifyEmbl(emblFile) {
  if (
    emblFile.search(/ID/) == -1 ||
    emblFile.search(/AC/) == -1 ||
    emblFile.search(/DE/) == -1 ||
    emblFile.search(/SQ/) == -1
  ) {
    alert("Please enter the contents of an EMBL file.");
    return false;
  }
  return true;
}

function verifyMaxDigits(theNumber, maxInput) {
  if (theNumber.search(/\d/) == -1) {
    alert("Please enter a number.");
    return false;
  }
  if (theNumber > maxInput) {
    alert("Please enter a number less than or equal to " + maxInput + ".");
    return false;
  }
}

function verifyDna(dnaSequence) {
  if (dnaSequence.search(/[^gatucryswkmbdhvnx\s]/i) != -1) {
    alert("The sequence contains non-DNA characters, which will be omitted.");
  }
  return true;
}

function verifyProtein(proteinSequence) {
  if (proteinSequence.search(/[^acdefghiklmnpqrstvwyz\*\s]/i) != -1) {
    alert(
      "The sequence contains non-protein characters, which will be omitted."
    );
  }
  return true;
}

function verifyGenBank(genBankFile) {
  if (
    genBankFile.search(/LOCUS/) == -1 ||
    genBankFile.search(/DEFINITION/) == -1 ||
    genBankFile.search(/ACCESSION/) == -1 ||
    genBankFile.search(/ORIGIN/) == -1
  ) {
    alert("Please enter the contents of a GenBank file.");
    return false;
  }
  return true;
}

function verifyGenBankFeat(genBankFile) {
  if (
    genBankFile.search(/LOCUS/) == -1 ||
    genBankFile.search(/DEFINITION/) == -1 ||
    genBankFile.search(/ACCESSION/) == -1 ||
    genBankFile.search(/ORIGIN/) == -1
  ) {
    alert("Please enter the contents of a GenBank file.");
    return false;
  }
  if (genBankFile.search(/FEATURES {13}/) == -1) {
    alert("The file has no defined features.");
    return false;
  }
  return true;
}

function verifyEmblFeat(emblFile) {
  if (
    emblFile.search(/ID/) == -1 ||
    emblFile.search(/AC/) == -1 ||
    emblFile.search(/DE/) == -1 ||
    emblFile.search(/SQ/) == -1
  ) {
    alert("Please enter the contents of an EMBL file.");
    return false;
  }
  if (emblFile.search(/^FT/m) == -1) {
    alert("The file has no defined features.");
    return false;
  }
  return true;
}

//deprecated
function writeGroupNum(
  text,
  tabIn,
  groupSize,
  basePerLine,
  startBase,
  stopBase
) {
  var i = parseInt(startBase);
  var k = 0;
  var lineOfText = "";
  var sepChar = " ";
  groupSize = parseInt(groupSize);
  basePerLine = parseInt(basePerLine);
  while (i < stopBase) {
    lineOfText = rightNum(i + 1, lineOfText, 8, tabIn);
    for (var j = 1; j <= basePerLine / groupSize; j++) {
      //makes a group each loop
      while (k < groupSize) {
        lineOfText = lineOfText + text.charAt(k + i);
        k = k + 1;
      }
      i = i + groupSize;
      k = 0;
      lineOfText = lineOfText + sepChar;
    }
    outputWindow.document.write(lineOfText + "\n");
    lineOfText = "";
  }
  return true;
}

function writeGroupNumDna(
  text,
  tabIn,
  groupSize,
  basePerLine,
  startBase,
  stopBase,
  strands,
  numberPosition
) {
  writeGroupNumDnaSetStart(
    text,
    tabIn,
    groupSize,
    basePerLine,
    startBase,
    stopBase,
    strands,
    numberPosition,
    0
  );
  return true;
}

function writeGroupNumDnaSetStart(
  text,
  tabIn,
  groupSize,
  basePerLine,
  startBase,
  stopBase,
  strands,
  numberPosition,
  numberingAdjustment
) {
  //strands can be 'one' or 'two'
  //numberPosition can be 'left', 'above', or 'right'
  var i = parseInt(startBase);
  var k = 0;
  var lineOfText = "";
  var lineNum = "";
  var sepChar = " ";
  var aboveNum = "";
  groupSize = parseInt(groupSize);
  basePerLine = parseInt(basePerLine);
  numberingAdjustment = parseInt(numberingAdjustment);
  while (i < stopBase) {
    lineNum = i + 1;
    for (var j = 1; j <= basePerLine / groupSize; j++) {
      //makes a group each loop
      while (k < groupSize) {
        if (i + k >= stopBase) {
          break;
        }
        lineOfText = lineOfText + text.charAt(k + i);
        k = k + 1;
      }
      lineOfText = lineOfText + sepChar;
      i = i + k;
      if (numberPosition == "above") {
        aboveNum =
          aboveNum +
          rightNum(
            adjustNumbering(i, numberingAdjustment),
            "",
            groupSize,
            tabIn
          );
      }
      if (i >= stopBase) {
        break;
      }
      k = 0;
    }
    if (numberPosition == "left") {
      outputWindow.document.write(
        rightNum(adjustNumbering(lineNum, numberingAdjustment), "", 8, tabIn) +
          lineOfText +
          "\n"
      );
      if (strands == "two") {
        outputWindow.document.write(
          rightNum(
            adjustNumbering(lineNum, numberingAdjustment),
            "",
            8,
            tabIn
          ) +
            complement(lineOfText) +
            "\n"
        );
        outputWindow.document.write("\n");
      }
    } else if (numberPosition == "right") {
      outputWindow.document.write(
        lineOfText + adjustNumbering(i, numberingAdjustment) + "\n"
      );
      if (strands == "two") {
        outputWindow.document.write(
          complement(lineOfText) +
            adjustNumbering(i, numberingAdjustment) +
            "\n"
        );
        outputWindow.document.write("\n");
      }
    } else if (numberPosition == "above") {
      outputWindow.document.write(aboveNum + "\n");
      outputWindow.document.write(lineOfText + "\n");
      if (strands == "two") {
        outputWindow.document.write(complement(lineOfText) + "\n");
        outputWindow.document.write("\n");
      }
    }
    aboveNum = "";
    lineOfText = "";
  }
  return true;

  //This function is used to change the numbering of base '0' to '1'
  function adjustNumbering(original, adjustment) {
    var adjusted = original + adjustment;
    if (adjustment < 0 && adjusted >= 0) {
      adjusted++;
    }
    return adjusted;
  }
}

function writeGroupNumProtein(
  text,
  tabIn,
  groupSize,
  basePerLine,
  startBase,
  stopBase,
  numberPosition
) {
  //numberPosition can be 'left', 'above', or 'right'
  var i = parseInt(startBase);
  var k = 0;
  var lineOfText = "";
  var lineNum = "";
  var sepChar = " ";
  var aboveNum = "";
  groupSize = parseInt(groupSize);
  basePerLine = parseInt(basePerLine);
  while (i < stopBase) {
    lineNum = i + 1;
    for (var j = 1; j <= basePerLine / groupSize; j++) {
      //makes a group each loop
      while (k < groupSize) {
        if (i + k >= stopBase) {
          break;
        }
        lineOfText = lineOfText + text.charAt(k + i);
        k = k + 1;
      }
      lineOfText = lineOfText + sepChar;
      i = i + k;
      if (numberPosition == "above") {
        aboveNum = aboveNum + rightNum(i, "", groupSize, tabIn);
      }
      if (i >= stopBase) {
        break;
      }
      k = 0;
    }
    if (numberPosition == "left") {
      outputWindow.document.write(
        rightNum(lineNum, "", 8, tabIn) + lineOfText + "\n"
      );
    } else if (numberPosition == "right") {
      outputWindow.document.write(lineOfText + i + "\n");
    } else if (numberPosition == "above") {
      outputWindow.document.write(aboveNum + "\n");
      outputWindow.document.write(lineOfText + "\n");
    }
    aboveNum = "";
    lineOfText = "";
  }
  return true;
}

function writeMutatedSequence(
  sequence,
  components,
  numMut,
  firstIndexToMutate,
  lastIndexToMutate
) {
  var currentChar = "";
  var randNum = 0;
  var maxNum = 0;
  var needNewChar = "";
  var componentsIndex = 0;

  numMut = parseInt(numMut);
  firstIndexToMutate = parseInt(firstIndexToMutate);
  lastIndexToMutate = parseInt(lastIndexToMutate);

  if (
    sequence.length <= firstIndexToMutate ||
    lastIndexToMutate < 0 ||
    lastIndexToMutate <= firstIndexToMutate
  ) {
    numMut = 0;
  }

  for (var i = 0; i < numMut; i++) {
    maxNum = sequence.length;
    randNum = Math.floor(Math.random() * maxNum);

    if (randNum < firstIndexToMutate || randNum > lastIndexToMutate) {
      numMut++;
      continue;
    }

    currentChar = sequence.charAt(randNum);
    needNewChar = true;
    while (needNewChar) {
      componentsIndex = Math.round(Math.random() * components.length);
      if (componentsIndex == components.length) {
        componentsIndex = 0;
      }
      if (components[componentsIndex] != currentChar) {
        needNewChar = false;
      }
    }
    sequence =
      sequence.substring(0, randNum) +
      components[componentsIndex] +
      sequence.substring(randNum + 1, sequence.length);
  }
  outputWindow.document.write(addReturns(sequence));
  return true;
}

function writeRandomSequence(components, lengthOut) {
  var sequence = "";
  var tempNum = 0;
  var tempChar = "";
  for (var j = 0; j < lengthOut; j++) {
    tempNum = Math.floor(Math.random() * components.length);
    tempChar = components[tempNum];
    sequence = sequence + tempChar;
    if (sequence.length == 60) {
      outputWindow.document.write(sequence + "\n");
      sequence = "";
    }
  }
  outputWindow.document.write(sequence + "\n");
  return true;
}

function writeRestrictionSites(sequence, arrayOfItems, dnaConformation) {
  var resultArray = new Array();
  var lookAhead = 50;
  var lowerLimit = 0;
  var upperLimit = sequence.length;
  var shiftValue = 0;
  var cutDistance;
  var matchExp;
  var matchPosition;
  var tempString;
  var backGroundClass;
  var matchArray;
  var timesFound = 0;
  if (dnaConformation == "circular") {
    shiftValue = sequence.substring(0, lookAhead).length;
    sequence =
      sequence.substring(sequence.length - lookAhead, sequence.length) +
      sequence +
      sequence.substring(0, lookAhead);
    lowerLimit = 0 + shiftValue;
    upperLimit = upperLimit + shiftValue;
  }
  outputWindow.document.write(
    '<table border="1" width="100%" cellspacing="0" cellpadding="2"><tbody>\n'
  );
  outputWindow.document.write(
    '<tr><td class="title" width="200px">' +
      "Site:" +
      '</td><td class="title">' +
      "Positions:" +
      "</td></tr>\n"
  );
  for (var i = 0; i < arrayOfItems.length; i++) {
    tempString = "none";
    backGroundClass = "many";
    matchExp = arrayOfItems[i].match(/\/.+\//) + "gi";
    matchPosition = 0;
    matchExp = eval(matchExp);
    cutDistance = parseFloat(
      arrayOfItems[i]
        .match(/\)\D*\d+/)
        .toString()
        .replace(/\)\D*/, "")
    );

    while ((matchArray = matchExp.exec(sequence))) {
      matchPosition = matchExp.lastIndex - cutDistance;
      if (matchPosition >= lowerLimit && matchPosition < upperLimit) {
        timesFound++;
        tempString = tempString + ", " + (matchPosition - shiftValue + 1);
      }
      matchExp.lastIndex = matchExp.lastIndex - RegExp.lastMatch.length + 1;
    }

    if (tempString.search(/\d/) != -1) {
      tempString = tempString.replace(/none,\s*/, "");
    }

    if (timesFound == 0) {
      backGroundClass = "none";
    } else if (timesFound == 1) {
      backGroundClass = "one";
    } else if (timesFound == 2) {
      backGroundClass = "two";
    } else if (timesFound == 3) {
      backGroundClass = "three";
    } else {
      backGroundClass = "many";
    }

    outputWindow.document.write(
      '<tr><td class="' +
        backGroundClass +
        '">' +
        arrayOfItems[i]
          .match(/\([^\(]+\)/)
          .toString()
          .replace(/\(|\)/g, "") +
        '</td><td class="' +
        backGroundClass +
        '">' +
        tempString +
        "</td></tr>\n"
    );

    timesFound = 0;
  }
  outputWindow.document.write("</tbody></table>\n");
  return true;
}

function writeSequenceStats(sequence, arrayOfItems) {
  //arrayOFItems are regular expressions. A number included with each regular expression serves as an adjustment for the percentage calculation. Any additional text will appear next to the pattern when the results are given.
  var originalLength = sequence.length;
  outputWindow.document.write(
    '<table border="1" width="100%" cellspacing="0" cellpadding="2"><tbody>\n'
  );
  outputWindow.document.write(
    '<tr><td class="title">' +
      "Pattern:" +
      '</td><td class="title">' +
      "Times found:" +
      '</td><td class="title">' +
      "Percentage:" +
      "</td></tr>\n"
  );
  for (var i = 0; i < arrayOfItems.length; i++) {
    var tempNumber = 0;
    var matchExp = arrayOfItems[i].match(/\/[^\/]+\//) + "gi";
    matchExp = eval(matchExp);
    if (sequence.search(matchExp) != -1) {
      tempNumber = sequence.match(matchExp).length;
    }
    var percentage = 0;
    if (originalLength + 1 - parseFloat(arrayOfItems[i].match(/\d+/)) > 0) {
      percentage =
        (100 * tempNumber) /
        (originalLength + 1 - parseFloat(arrayOfItems[i].match(/\d+/)));
    }

    outputWindow.document.write(
      "<tr><td>" +
        arrayOfItems[i]
          .match(/\([^\(]+\)\b/)
          .toString()
          .replace(/\(|\)/g, "") +
        "</td><td>" +
        tempNumber +
        "</td><td>" +
        percentage.toFixed(2) +
        "</td></tr>\n"
    );
  }
  outputWindow.document.write("</tbody></table>\n");
  return true;
}

function writeShuffledSequence(sequence) {
  var tempSeq = "";
  var tempChar = "";
  var tempString1 = "";
  var tempString2 = "";
  var randNum = 0;
  var maxNum = 0;
  while (sequence.length > 0) {
    maxNum = sequence.length;
    randNum = Math.floor(Math.random() * maxNum);
    tempChar = sequence.charAt(randNum);
    tempSeq = tempSeq + tempChar;
    tempString1 = sequence.substring(0, randNum);
    tempString2 = sequence.substring(randNum + 1, sequence.length);
    sequence = tempString1 + tempString2;
    if (tempSeq.length == 60) {
      outputWindow.document.write(tempSeq + "\n");
      tempSeq = "";
    }
  }
  outputWindow.document.write(tempSeq + "\n");
  return true;
}
