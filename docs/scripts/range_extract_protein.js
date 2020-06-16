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

function rangeExtract(theDocument) {
  var newDna = "";
  var maxInput = 500000000;
  var matchFound = false;
  var ranges = new Array();

  if (testScript() == false) {
    return false;
  }

  if (
    checkFormElement(theDocument.forms[0].elements[0]) == false ||
    checkSequenceLength(theDocument.forms[0].elements[0].value, maxInput) ==
      false ||
    checkFormElement(theDocument.forms[0].elements[1]) == false
  ) {
    return false;
  }

  var arrayOfRanges = theDocument.forms[0].elements[1].value.split(/,/);
  var arrayOfStartAndEnd;
  for (var i = 0; i < arrayOfRanges.length; i++) {
    arrayOfStartAndEnd = arrayOfRanges[i].split(/\.\./);
    if (arrayOfStartAndEnd.length == 1) {
      matchFound = true;
      ranges.push(new Range(arrayOfStartAndEnd[0], arrayOfStartAndEnd[0]));
    } else if (arrayOfStartAndEnd.length == 2) {
      matchFound = true;
      ranges.push(new Range(arrayOfStartAndEnd[0], arrayOfStartAndEnd[1]));
    }
  }
  if (matchFound == false) {
    alert("No ranges were entered.");
    return false;
  }

  openWindow("Range Extractor Protein");
  openPre();
  var arrayOfFasta = getArrayOfFasta(theDocument.forms[0].elements[0].value);
  for (var i = 0; i < arrayOfFasta.length; i++) {
    newProtein = getSequenceFromFasta(arrayOfFasta[i]);
    title = getTitleFromFasta(arrayOfFasta[i]);
    newProtein = removeNonProteinAllowDegen(newProtein);
    outputWindow.document.write(
      getFastaTitleFromTitleAndSequence(title, newProtein)
    );
    writeSequenceRanges(
      newProtein,
      ranges,
      theDocument.forms[0].elements[5].options[
        theDocument.forms[0].elements[5].selectedIndex
      ].value
    );
  }
  closePre();
  closeWindow();
  return true;
}

function writeSequenceRanges(sequence, ranges, segmentType) {
  var rangeGroup = new RangeGroup(segmentType);

  //replace values like 'start' and 'end' and 'length' with numerical values
  //replace values like '(end - 3)' with numerical values
  var center_base = Math.round(sequence.length / 2);

  for (var i = 0; i < ranges.length; i++) {
    ranges[i].start = ranges[i].start.toString().replace(/start|begin/gi, 1);
    ranges[i].start = ranges[i].start
      .toString()
      .replace(/stop|end/gi, sequence.length);
    ranges[i].start = ranges[i].start
      .toString()
      .replace(/length/gi, sequence.length);
    ranges[i].start = ranges[i].start
      .toString()
      .replace(/middle|center|centre/gi, center_base);

    ranges[i].stop = ranges[i].stop.toString().replace(/start|begin/gi, 1);
    ranges[i].stop = ranges[i].stop
      .toString()
      .replace(/stop|end/gi, sequence.length);
    ranges[i].stop = ranges[i].stop
      .toString()
      .replace(/length/gi, sequence.length);
    ranges[i].stop = ranges[i].stop
      .toString()
      .replace(/middle|center|centre/gi, center_base);

    try {
      ranges[i].start = parseInt(eval(ranges[i].start.toString()));
    } catch (e) {
      alert("Could not evaluate the following expression: " + ranges[i].start);
      return false;
    }
    try {
      ranges[i].stop = parseInt(eval(ranges[i].stop.toString()));
    } catch (e) {
      alert("Could not evaluate the following expression: " + ranges[i].stop);
      return false;
    }

    rangeGroup.addRange(ranges[i]);
  }
  rangeGroup.writeRanges(sequence);
}

//Range class getSequence() method
function getSequence(sequence) {
  var problem = false;
  var warnings = new Array();

  if (this.start == this.stop) {
    if (this.start < 1) {
      problem = true;
      warnings.push("position value of " + this.start + " is less than 1");
    }

    if (this.start > sequence.length) {
      problem = true;
      warnings.push(
        "position value of " +
          this.start +
          " is greater than the sequence length"
      );
    }
  } else {
    if (this.start < 1) {
      problem = true;
      warnings.push("position value of " + this.start + " is less than 1");
    }
    if (this.stop < 1) {
      problem = true;
      warnings.push("position value of " + this.stop + " is less than 1");
    }
    if (this.start > sequence.length) {
      problem = true;
      warnings.push(
        "position value of " +
          this.start +
          " is greater than the sequence length"
      );
    }
    if (this.stop > sequence.length) {
      problem = true;
      warnings.push(
        "position value of " +
          this.stop +
          " is greater than the sequence length"
      );
    }
    if (this.start > this.stop) {
      problem = true;
      warnings.push(
        "stop position is less than start position in range " +
          this.start +
          " to " +
          this.stop
      );
    }
  }

  if (problem) {
    alert(
      "An entry was skipped because of the following:\n" + warnings.join(",\n")
    );
    return false;
  } else if (this.start == this.stop) {
    return sequence.charAt(this.start - 1);
  } else {
    return sequence.substring(this.start - 1, this.stop);
  }
}

//Range class getTitle() method
function getTitle() {
  if (this.start == this.stop) {
    return "&gt;residue " + this.start;
  } else {
    return "&gt;residue " + this.start + " to " + this.stop;
  }
}

//Range class
function Range(start, stop) {
  this.start = start;
  this.stop = stop;
}

//create and throw away a prototype object
new Range(0, 0);

// define object methods
Range.prototype.getSequence = getSequence;
Range.prototype.getTitle = getTitle;

//RangeGroup class addRange() method
function addRange(range) {
  this.ranges.push(range);
}

//RangeGroup class writeRanges() method
function writeRanges(sequence) {
  var sequenceArray = new Array();

  if (this.type == "new_sequence") {
    for (var i = 0; i < this.ranges.length; i++) {
      if (this.ranges[i].getSequence(sequence) != "") {
        sequenceArray.push(this.ranges[i].getSequence(sequence));
      }
    }
    outputWindow.document.write(addReturns(sequenceArray.join("")) + "\n\n");
    return true;
  }

  if (this.type == "separate") {
    for (var i = 0; i < this.ranges.length; i++) {
      if (this.ranges[i].getSequence(sequence) != "") {
        outputWindow.document.write(this.ranges[i].getTitle() + "\n");
        outputWindow.document.write(
          addReturns(this.ranges[i].getSequence(sequence)) + "\n\n"
        );
      }
    }
    return true;
  }

  if (this.type == "uppercased") {
    var re;
    sequence = sequence.toLowerCase();
    for (var i = 0; i < this.ranges.length; i++) {
      if (this.ranges[i].getSequence(sequence) != "") {
        if (this.ranges[i].start > 1) {
          re =
            "(.{" +
            (this.ranges[i].start - 1) +
            "})\\B(.{" +
            (this.ranges[i].stop - this.ranges[i].start + 1) +
            "})";
          re = new RegExp(re);
          sequence = sequence.replace(re, function (str, p1, p2, offset, s) {
            return p1 + p2.toUpperCase();
          });
        } else {
          re = "(.{" + (this.ranges[i].stop - this.ranges[i].start + 1) + "})";
          re = new RegExp(re);
          sequence = sequence.replace(re, function (str, p1, offset, s) {
            return p1.toUpperCase();
          });
        }
      }
    }
    outputWindow.document.write(addReturns(sequence) + "\n\n");
    return true;
  }

  if (this.type == "lowercased") {
    var re;
    sequence = sequence.toUpperCase();
    for (var i = 0; i < this.ranges.length; i++) {
      if (this.ranges[i].getSequence(sequence) != "") {
        if (this.ranges[i].start > 1) {
          re =
            "(.{" +
            (this.ranges[i].start - 1) +
            "})\\B(.{" +
            (this.ranges[i].stop - this.ranges[i].start + 1) +
            "})";
          re = new RegExp(re);
          sequence = sequence.replace(re, function (str, p1, p2, offset, s) {
            return p1 + p2.toLowerCase();
          });
        } else {
          re = "(.{" + (this.ranges[i].stop - this.ranges[i].start + 1) + "})";
          re = new RegExp(re);
          sequence = sequence.replace(re, function (str, p1, offset, s) {
            return p1.toLowerCase();
          });
        }
      }
    }
    outputWindow.document.write(addReturns(sequence) + "\n\n");
    return true;
  }

  if (this.type == "randomized") {
    var re;
    for (var i = 0; i < this.ranges.length; i++) {
      if (this.ranges[i].getSequence(sequence) != "") {
        if (this.ranges[i].start > 1) {
          re =
            "(.{" +
            (this.ranges[i].start - 1) +
            "})\\B(.{" +
            (this.ranges[i].stop - this.ranges[i].start + 1) +
            "})";
          re = new RegExp(re);
          sequence = sequence.replace(re, function (str, p1, p2, offset, s) {
            return (
              p1 +
              getRandomSequence(
                [
                  "A",
                  "C",
                  "D",
                  "E",
                  "F",
                  "G",
                  "H",
                  "I",
                  "K",
                  "L",
                  "M",
                  "N",
                  "P",
                  "Q",
                  "R",
                  "S",
                  "T",
                  "V",
                  "W",
                  "Y",
                ],
                p2.length
              )
            );
          });
        } else {
          re = "(.{" + (this.ranges[i].stop - this.ranges[i].start + 1) + "})";
          re = new RegExp(re);
          sequence = sequence.replace(re, function (str, p1, offset, s) {
            return getRandomSequence(
              [
                "A",
                "C",
                "D",
                "E",
                "F",
                "G",
                "H",
                "I",
                "K",
                "L",
                "M",
                "N",
                "P",
                "Q",
                "R",
                "S",
                "T",
                "V",
                "W",
                "Y",
              ],
              p1.length
            );
          });
        }
      }
    }
    outputWindow.document.write(addReturns(sequence) + "\n\n");
    return true;
  }

  if (this.type == "preserved") {
    var re;
    var randomSequence = getRandomSequence(
      [
        "A",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y",
      ],
      sequence.length
    );
    for (var i = 0; i < this.ranges.length; i++) {
      if (this.ranges[i].getSequence(sequence) != "") {
        if (this.ranges[i].start > 1) {
          re =
            "(.{" +
            (this.ranges[i].start - 1) +
            "})\\B(.{" +
            (this.ranges[i].stop - this.ranges[i].start + 1) +
            "})";
          re = new RegExp(re);
          randomSequence = randomSequence.replace(re, function (
            str,
            p1,
            p2,
            offset,
            s
          ) {
            return p1 + sequence.substring(p1.length, p1.length + p2.length);
          });
        } else {
          re = "(.{" + (this.ranges[i].stop - this.ranges[i].start + 1) + "})";
          re = new RegExp(re);
          randomSequence = randomSequence.replace(re, function (
            str,
            p1,
            offset,
            s
          ) {
            return sequence.substring(offset, offset + p1.length);
          });
        }
      }
    }
    outputWindow.document.write(addReturns(randomSequence) + "\n\n");
    return true;
  }
}

//RangeGroup class
function RangeGroup(type) {
  this.type = type;
  this.ranges = new Array();
}

//create and throw away a prototype object
new RangeGroup("", "");

// define object methods
RangeGroup.prototype.addRange = addRange;
RangeGroup.prototype.writeRanges = writeRanges;
