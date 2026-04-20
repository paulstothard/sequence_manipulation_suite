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

// Written by Paul Stothard, University of Alberta, Canada
// Codon-level affine-gap alignment using the Schneider et al. (2005) matrix.

(function (global) {
  "use strict";

  var NEG_INF = -1073741824;
  var SMALL_PROBLEM_SIZE = 4096;

  var STATE_M = 0;
  var STATE_X = 1;
  var STATE_Y = 2;
  var MODE_START = 3;
  var MODE_ANY = 4;

  var CODON_ORDER = [
    "AAA","AAC","AAG","AAT","ACA","ACC","ACG","ACT",
    "AGA","AGC","AGG","AGT","ATA","ATC","ATG","ATT",
    "CAA","CAC","CAG","CAT","CCA","CCC","CCG","CCT",
    "CGA","CGC","CGG","CGT","CTA","CTC","CTG","CTT",
    "GAA","GAC","GAG","GAT","GCA","GCC","GCG","GCT",
    "GGA","GGC","GGG","GGT","GTA","GTC","GTG","GTT",
    "TAA","TAC","TAG","TAT","TCA","TCC","TCG","TCT",
    "TGA","TGC","TGG","TGT","TTA","TTC","TTG","TTT"
  ];

  var CODON_LOOKUP = (function () {
    var tbl = {};
    for (var i = 0; i < CODON_ORDER.length; i++) {
      tbl[CODON_ORDER[i]] = i;
    }
    return tbl;
  })();

  var CODON_SCORES = new Int16Array([
    11, -2, 9, -1, -2, -6, -3, -5, 5, -5, 3, -4, -6, -13, -7, -11, 0, -6, -1, -5, -8, -11, -8, -10, 2, 0, 1, 0, -10, -13, -13, -12, -2, -8, -5, -8, -6, -9, -7, -9, -7, -10, -8, -9, -8, -12, -11, -11, -50, -14, -50, -13, -7, -10, -8, -9, -50, -13, -13, -12, -10, -18, -11, -17,
    -2, 13, -3, 10, -3, 0, -3, -1, -5, 4, -5, 3, -10, -7, -9, -9, -5, 0, -5, -1, -10, -8, -9, -9, -8, -5, -7, -6, -13, -11, -14, -12, -6, 0, -6, -1, -7, -5, -6, -6, -5, -1, -4, -3, -10, -8, -11, -10, -50, -6, -50, -7, -6, -4, -6, -5, -50, -7, -16, -8, -13, -12, -13, -13,
    9, -3, 11, -2, -4, -6, -3, -6, 3, -5, 4, -5, -8, -13, -5, -12, -1, -6, 0, -6, -10, -11, -8, -11, 1, 0, 2, 0, -11, -14, -11, -13, -4, -9, -3, -8, -8, -10, -7, -10, -9, -10, -7, -9, -10, -13, -10, -12, -50, -14, -50, -14, -8, -10, -8, -10, -50, -13, -11, -13, -12, -19, -11, -17,
    -1, 10, -2, 12, -2, -2, -2, 0, -4, 2, -5, 5, -9, -9, -9, -7, -4, -1, -5, 1, -9, -8, -8, -7, -7, -6, -7, -4, -12, -12, -13, -10, -4, 0, -5, 1, -6, -6, -6, -4, -5, -3, -5, -1, -9, -9, -11, -8, -50, -7, -50, -5, -5, -5, -5, -4, -50, -8, -15, -5, -11, -14, -12, -11,
    -2, -3, -4, -2, 11, 9, 10, 9, -3, 0, -5, 0, 0, -4, 0, -3, -4, -8, -6, -7, -1, -4, -2, -3, -8, -9, -8, -8, -6, -9, -7, -8, -6, -10, -7, -8, 2, 0, 0, 0, -5, -7, -6, -6, 0, -3, -2, -3, -50, -14, -50, -12, 2, 0, 1, 1, -50, -9, -13, -7, -4, -12, -6, -11,
    -6, 0, -6, -2, 9, 12, 9, 9, -7, 2, -7, 0, -2, -1, -3, -3, -7, -6, -8, -7, -4, -1, -4, -3, -9, -7, -9, -8, -9, -7, -9, -9, -9, -7, -9, -9, 0, 2, 0, 0, -8, -4, -7, -6, -3, -1, -4, -3, -50, -11, -50, -11, 0, 2, 0, 0, -50, -6, -15, -7, -8, -9, -8, -11,
    -3, -3, -3, -2, 10, 9, 12, 9, -4, 0, -3, 0, -1, -4, 0, -3, -5, -8, -5, -8, -3, -3, -2, -4, -7, -7, -6, -8, -7, -8, -6, -7, -7, -9, -7, -9, 0, 0, 2, 0, -6, -6, -5, -6, -1, -3, -1, -2, -50, -13, -50, -11, 1, 0, 2, 0, -50, -8, -12, -7, -6, -11, -5, -10,
    -5, -1, -6, 0, 9, 9, 9, 11, -6, 0, -7, 2, -2, -3, -2, 0, -6, -8, -8, -6, -3, -3, -3, -1, -9, -8, -9, -7, -8, -8, -8, -6, -8, -8, -8, -6, 0, 0, 0, 2, -7, -6, -6, -4, -2, -2, -3, 0, -50, -12, -50, -10, 1, 0, 1, 2, -50, -7, -16, -5, -7, -10, -7, -9,
    5, -5, 3, -4, -3, -7, -4, -6, 13, -2, 11, -1, -5, -12, -7, -11, 0, -4, -3, -4, -10, -11, -8, -11, 10, 7, 9, 8, -9, -12, -11, -11, -6, -11, -8, -11, -7, -10, -8, -9, -1, -6, -4, -6, -7, -12, -11, -11, -50, -14, -50, -13, -8, -10, -9, -9, -50, -9, -7, -8, -10, -17, -11, -16,
    -5, 4, -5, 2, 0, 2, 0, 0, -2, 12, -2, 11, -8, -6, -8, -7, -6, -2, -6, -4, -8, -5, -7, -7, -5, -2, -5, -4, -12, -10, -12, -11, -7, -3, -7, -4, -4, -1, -3, -3, -1, 3, 0, 0, -8, -6, -9, -7, -50, -9, -50, -9, -2, 0, -1, -1, -50, 0, -13, -2, -11, -11, -12, -12,
    3, -5, 4, -5, -5, -7, -3, -7, 11, -2, 13, -2, -7, -12, -5, -11, -2, -4, -1, -4, -10, -11, -8, -12, 9, 8, 10, 8, -10, -11, -10, -12, -8, -11, -6, -11, -8, -9, -7, -10, -4, -6, -1, -7, -9, -11, -9, -12, -50, -14, -50, -12, -10, -10, -8, -10, -50, -9, -4, -9, -11, -18, -11, -16,
    -4, 3, -5, 5, 0, 0, 0, 2, -1, 11, -2, 13, -7, -8, -8, -5, -5, -4, -6, -2, -7, -7, -7, -6, -5, -5, -5, -2, -12, -12, -13, -10, -6, -4, -6, -2, -3, -3, -3, -2, -1, 0, 0, 3, -8, -8, -9, -5, -50, -10, -50, -7, -1, -2, -1, 0, -50, -2, -12, 0, -11, -13, -11, -11,
    -6, -10, -8, -9, 0, -2, -1, -2, -5, -8, -7, -7, 13, 9, 3, 9, -8, -12, -10, -10, -7, -10, -8, -9, -9, -11, -10, -10, 2, 0, 0, 0, -9, -15, -11, -13, -3, -6, -4, -5, -9, -11, -10, -11, 6, 3, 3, 3, -50, -13, -50, -11, -5, -8, -7, -8, -50, -12, -14, -11, 2, -6, 0, -5,
    -13, -7, -13, -9, -4, -1, -4, -3, -12, -6, -12, -8, 9, 12, 0, 10, -12, -11, -13, -12, -12, -10, -11, -12, -15, -13, -14, -13, -1, 1, -1, 0, -14, -14, -15, -16, -7, -4, -6, -6, -14, -10, -12, -13, 2, 6, 2, 3, -50, -11, -50, -11, -10, -9, -10, -10, -50, -10, -16, -11, -1, -3, -2, -5,
    -7, -9, -5, -9, 0, -3, 0, -2, -7, -8, -5, -8, 3, 0, 14, 1, -7, -11, -6, -9, -8, -10, -7, -9, -10, -10, -8, -10, 1, 0, 1, 0, -11, -14, -9, -14, -3, -5, -3, -5, -10, -11, -8, -11, 0, -1, 1, -1, -50, -12, -50, -11, -5, -8, -4, -7, -50, -12, -10, -11, 0, -6, 2, -6,
    -11, -9, -12, -7, -3, -3, -3, 0, -11, -7, -11, -5, 9, 10, 1, 12, -12, -12, -12, -9, -11, -11, -11, -9, -12, -14, -14, -11, -1, 0, -1, 1, -13, -15, -13, -12, -6, -6, -6, -4, -13, -12, -13, -10, 2, 3, 2, 5, -50, -12, -50, -10, -9, -9, -9, -8, -50, -11, -15, -9, -1, -5, -1, -3,
    0, -5, -1, -4, -4, -7, -5, -6, 0, -6, -2, -5, -8, -12, -7, -12, 12, 2, 10, 3, 0, -3, 0, -3, 2, 0, 0, 0, -3, -7, -5, -6, 0, -6, -1, -6, -5, -8, -5, -7, -6, -9, -7, -9, -7, -10, -9, -9, -50, -8, -50, -7, -4, -6, -5, -6, -50, -9, -9, -8, -6, -12, -6, -11,
    -6, 0, -6, -1, -8, -6, -8, -8, -4, -2, -4, -4, -12, -11, -11, -12, 2, 14, 1, 12, -5, -2, -4, -4, -1, 3, -1, 1, -8, -4, -8, -6, -8, -4, -7, -6, -10, -8, -8, -9, -10, -7, -10, -9, -11, -10, -12, -11, -50, 2, -50, 0, -7, -5, -7, -6, -50, -4, -11, -5, -9, -5, -9, -7,
    -1, -5, 0, -5, -6, -8, -5, -8, -3, -6, -1, -6, -10, -13, -6, -12, 10, 1, 11, 2, -2, -4, 0, -4, 0, 0, 2, 0, -5, -6, -4, -6, -2, -7, 0, -7, -7, -8, -5, -8, -9, -9, -7, -9, -9, -11, -9, -10, -50, -8, -50, -8, -6, -7, -5, -7, -50, -10, -7, -10, -7, -13, -6, -13,
    -5, -1, -6, 1, -7, -7, -8, -6, -4, -4, -4, -2, -10, -12, -9, -9, 3, 12, 2, 14, -4, -3, -4, -2, -1, 0, -1, 3, -7, -5, -8, -3, -6, -5, -7, -4, -8, -10, -8, -8, -9, -8, -9, -6, -11, -11, -11, -10, -50, 0, -50, 2, -7, -7, -7, -5, -50, -5, -9, -3, -8, -7, -8, -5,
    -8, -10, -10, -9, -1, -4, -3, -3, -10, -8, -10, -7, -7, -12, -8, -11, 0, -5, -2, -4, 12, 10, 11, 10, -6, -8, -6, -7, -2, -7, -5, -5, -8, -12, -9, -11, 0, -3, -2, -3, -9, -9, -9, -9, -5, -9, -8, -8, -50, -15, -50, -14, 2, -1, 0, 0, -50, -13, -14, -11, -5, -13, -6, -12,
    -11, -8, -11, -8, -4, -1, -3, -3, -11, -5, -11, -7, -10, -10, -10, -11, -3, -2, -4, -3, 10, 13, 10, 10, -8, -5, -8, -6, -6, -3, -7, -5, -11, -10, -10, -11, -3, -1, -2, -2, -10, -8, -9, -9, -8, -6, -9, -8, -50, -11, -50, -12, 0, 2, 0, 0, -50, -10, -17, -11, -8, -9, -8, -11,
    -8, -9, -8, -8, -2, -4, -2, -3, -8, -7, -8, -7, -8, -11, -7, -11, 0, -4, 0, -4, 11, 10, 13, 10, -5, -5, -3, -6, -3, -6, -3, -5, -9, -11, -8, -11, -1, -2, 0, -3, -9, -8, -7, -9, -7, -8, -7, -8, -50, -13, -50, -13, 0, 0, 1, 0, -50, -11, -11, -10, -6, -12, -5, -11,
    -10, -9, -11, -7, -3, -3, -4, -1, -11, -7, -12, -6, -9, -12, -9, -9, -3, -4, -4, -2, 10, 10, 10, 12, -8, -7, -9, -4, -6, -6, -7, -2, -10, -12, -11, -10, -2, -2, -2, 0, -10, -9, -10, -8, -8, -9, -9, -6, -50, -13, -50, -11, 0, 0, 0, 2, -50, -10, -17, -8, -8, -12, -7, -9,
    2, -8, 1, -7, -8, -9, -7, -9, 10, -5, 9, -5, -9, -15, -10, -12, 2, -1, 0, -1, -6, -8, -5, -8, 13, 11, 11, 12, -6, -9, -9, -9, -8, -13, -10, -13, -9, -11, -9, -12, -5, -9, -6, -8, -10, -14, -13, -13, -50, -11, -50, -9, -8, -10, -9, -11, -50, -6, -5, -6, -10, -16, -10, -15,
    0, -5, 0, -6, -9, -7, -7, -8, 7, -2, 8, -5, -11, -13, -10, -14, 0, 3, 0, 0, -8, -5, -5, -7, 11, 15, 11, 12, -9, -5, -9, -7, -11, -10, -10, -12, -11, -8, -8, -11, -9, -5, -8, -8, -13, -11, -12, -13, -50, -7, -50, -8, -11, -7, -8, -10, -50, -1, -7, -4, -11, -11, -11, -15,
    1, -7, 2, -7, -8, -9, -6, -9, 9, -5, 10, -5, -10, -14, -8, -14, 0, -1, 2, -1, -6, -8, -3, -9, 11, 11, 13, 11, -7, -8, -6, -8, -10, -12, -7, -12, -9, -11, -7, -10, -7, -8, -4, -9, -11, -12, -10, -12, -50, -11, -50, -11, -9, -9, -7, -10, -50, -6, -2, -6, -9, -15, -8, -14,
    0, -6, 0, -4, -8, -8, -8, -7, 8, -4, 8, -2, -10, -13, -10, -11, 0, 1, 0, 3, -7, -6, -6, -4, 12, 12, 11, 14, -8, -7, -9, -5, -9, -11, -10, -10, -9, -10, -9, -9, -8, -7, -8, -5, -12, -12, -12, -10, -50, -7, -50, -5, -9, -8, -8, -7, -50, -3, -7, -1, -10, -12, -9, -11,
    -10, -13, -11, -12, -6, -9, -7, -8, -9, -12, -10, -12, 2, -1, 1, -1, -3, -8, -5, -7, -2, -6, -3, -6, -6, -9, -7, -8, 11, 7, 8, 8, -12, -17, -13, -16, -6, -9, -7, -8, -12, -15, -13, -14, 0, -3, -1, -3, -50, -10, -50, -9, -4, -8, -5, -8, -50, -12, -9, -10, 9, -3, 8, -2,
    -13, -11, -14, -12, -9, -7, -8, -8, -12, -10, -11, -12, 0, 1, 0, 0, -7, -4, -6, -5, -7, -3, -6, -6, -9, -5, -8, -7, 7, 11, 7, 9, -15, -15, -14, -17, -9, -7, -8, -9, -15, -12, -14, -14, -2, 0, -2, -2, -50, -7, -50, -8, -9, -7, -9, -9, -50, -8, -11, -9, 6, 0, 6, -1,
    -13, -14, -11, -13, -7, -9, -6, -8, -11, -12, -10, -13, 0, -1, 1, -1, -5, -8, -4, -8, -5, -7, -3, -7, -9, -9, -6, -9, 8, 7, 10, 7, -14, -17, -13, -17, -8, -9, -6, -9, -15, -14, -12, -14, -1, -3, 0, -3, -50, -10, -50, -9, -7, -9, -6, -9, -50, -12, -8, -11, 7, -3, 8, -2,
    -12, -12, -13, -10, -8, -9, -7, -6, -11, -11, -12, -10, 0, 0, 0, 1, -6, -6, -6, -3, -5, -5, -5, -2, -9, -7, -8, -5, 8, 9, 7, 11, -14, -15, -14, -14, -8, -8, -7, -7, -14, -13, -13, -12, -2, -2, -2, 0, -50, -8, -50, -7, -8, -9, -8, -6, -50, -9, -11, -8, 6, -1, 6, 0,
    -2, -6, -4, -4, -6, -9, -7, -8, -6, -7, -8, -6, -9, -14, -11, -13, 0, -8, -2, -6, -8, -11, -9, -10, -8, -11, -10, -9, -12, -15, -14, -14, 11, 2, 9, 3, -3, -7, -4, -6, -2, -6, -3, -5, -6, -10, -8, -9, -50, -15, -50, -13, -7, -10, -8, -9, -50, -16, -17, -14, -12, -18, -12, -17,
    -8, 0, -9, 0, -10, -7, -9, -8, -11, -3, -11, -4, -15, -14, -14, -15, -6, -4, -7, -5, -12, -10, -11, -12, -13, -10, -12, -11, -17, -15, -17, -15, 2, 12, 3, 10, -7, -5, -6, -7, -5, -1, -5, -3, -11, -9, -13, -10, -50, -10, -50, -11, -9, -8, -9, -9, -50, -12, -20, -13, -16, -16, -17, -18,
    -5, -6, -3, -5, -7, -9, -7, -8, -8, -7, -6, -6, -11, -15, -9, -13, -1, -7, 0, -7, -9, -10, -8, -11, -10, -10, -7, -10, -13, -14, -13, -14, 9, 3, 10, 3, -4, -6, -2, -6, -5, -6, -2, -6, -7, -10, -7, -10, -50, -15, -50, -14, -8, -10, -8, -10, -50, -15, -15, -15, -13, -18, -12, -18,
    -8, -1, -8, 1, -8, -9, -9, -6, -11, -4, -11, -2, -13, -16, -14, -12, -6, -6, -7, -4, -11, -11, -11, -10, -13, -12, -12, -10, -16, -17, -17, -14, 3, 10, 3, 12, -6, -7, -6, -4, -5, -4, -5, -1, -10, -11, -12, -8, -50, -12, -50, -8, -9, -10, -9, -7, -50, -14, -19, -11, -16, -19, -15, -15,
    -6, -7, -8, -6, 2, 0, 0, 0, -7, -4, -8, -3, -3, -7, -3, -6, -5, -10, -7, -8, 0, -3, -1, -2, -9, -11, -9, -9, -6, -9, -8, -8, -3, -7, -4, -6, 11, 8, 9, 9, -1, -3, -1, -2, 1, -2, 0, -1, -50, -14, -50, -12, 2, 0, 1, 0, -50, -8, -14, -7, -5, -12, -6, -10,
    -9, -5, -10, -6, 0, 2, 0, 0, -10, -1, -9, -3, -6, -4, -5, -6, -8, -8, -8, -10, -3, -1, -2, -2, -11, -8, -11, -10, -9, -7, -9, -8, -7, -5, -6, -7, 8, 11, 8, 9, -3, 0, -3, -3, -2, 1, -2, -1, -50, -12, -50, -12, 0, 2, 0, 0, -50, -6, -14, -7, -8, -9, -8, -11,
    -7, -6, -7, -6, 0, 0, 2, 0, -8, -3, -7, -3, -4, -6, -3, -6, -5, -8, -5, -8, -2, -2, 0, -2, -9, -8, -7, -9, -7, -8, -6, -7, -4, -6, -2, -6, 9, 8, 12, 8, -2, -2, 0, -2, 0, -1, 1, 0, -50, -11, -50, -12, 1, 0, 3, 0, -50, -7, -12, -7, -7, -10, -5, -10,
    -9, -6, -10, -4, 0, 0, 0, 2, -9, -3, -10, -2, -5, -6, -5, -4, -7, -9, -8, -8, -3, -2, -3, 0, -12, -11, -10, -9, -8, -9, -9, -7, -6, -7, -6, -4, 9, 9, 8, 11, -3, -3, -3, -1, -1, -1, -2, 1, -50, -12, -50, -10, 0, 0, 0, 2, -50, -7, -16, -5, -7, -11, -7, -9,
    -7, -5, -9, -5, -5, -8, -6, -7, -1, -1, -4, -1, -9, -14, -10, -13, -6, -10, -9, -9, -9, -10, -9, -10, -5, -9, -7, -8, -12, -15, -15, -14, -2, -5, -5, -5, -1, -3, -2, -3, 12, 9, 11, 10, -4, -9, -8, -8, -50, -18, -50, -15, -5, -7, -6, -7, -50, -9, -11, -8, -12, -17, -13, -15,
    -10, -1, -10, -3, -7, -4, -6, -6, -6, 3, -6, 0, -11, -10, -11, -12, -9, -7, -9, -8, -9, -8, -8, -9, -9, -5, -8, -7, -15, -12, -14, -13, -6, -1, -6, -4, -3, 0, -2, -3, 9, 12, 9, 10, -8, -5, -9, -7, -50, -12, -50, -13, -6, -5, -5, -6, -50, -4, -12, -5, -14, -14, -13, -14,
    -8, -4, -7, -5, -6, -7, -5, -6, -4, 0, -1, 0, -10, -12, -8, -13, -7, -10, -7, -9, -9, -9, -7, -10, -6, -8, -4, -8, -13, -14, -12, -13, -3, -5, -2, -5, -1, -3, 0, -3, 11, 9, 12, 9, -5, -8, -5, -7, -50, -16, -50, -14, -5, -7, -4, -7, -50, -8, -6, -8, -12, -16, -10, -15,
    -9, -3, -9, -1, -6, -6, -6, -4, -6, 0, -7, 3, -11, -13, -11, -10, -9, -9, -9, -6, -9, -9, -9, -8, -8, -8, -9, -5, -14, -14, -14, -12, -5, -3, -6, -1, -2, -3, -2, -1, 10, 10, 9, 13, -7, -8, -8, -5, -50, -14, -50, -11, -6, -6, -5, -5, -50, -6, -13, -3, -13, -16, -13, -13,
    -8, -10, -10, -9, 0, -3, -1, -2, -7, -8, -9, -8, 6, 2, 0, 2, -7, -11, -9, -11, -5, -8, -7, -8, -10, -13, -11, -12, 0, -2, -1, -2, -6, -11, -7, -10, 1, -2, 0, -1, -4, -8, -5, -7, 11, 8, 10, 9, -50, -14, -50, -12, -4, -7, -5, -6, -50, -11, -14, -9, 1, -7, 0, -6,
    -12, -8, -13, -9, -3, -1, -3, -2, -12, -6, -11, -8, 3, 6, -1, 3, -10, -10, -11, -11, -9, -6, -8, -9, -14, -11, -12, -12, -3, 0, -3, -2, -10, -9, -10, -11, -2, 1, -1, -1, -9, -5, -8, -8, 8, 12, 8, 9, -50, -10, -50, -11, -7, -5, -7, -7, -50, -8, -15, -9, -3, -3, -3, -5,
    -11, -11, -10, -11, -2, -4, -1, -3, -11, -9, -9, -9, 3, 2, 1, 2, -9, -12, -9, -11, -8, -9, -7, -9, -13, -12, -10, -12, -1, -2, 0, -2, -8, -13, -7, -12, 0, -2, 1, -2, -8, -9, -5, -8, 10, 8, 11, 8, -50, -13, -50, -13, -6, -8, -5, -7, -50, -11, -12, -10, -1, -7, 0, -7,
    -11, -10, -12, -8, -3, -3, -2, 0, -11, -7, -12, -5, 3, 3, -1, 5, -9, -11, -10, -10, -8, -8, -8, -6, -13, -13, -12, -10, -3, -2, -3, 0, -9, -10, -10, -8, -1, -1, 0, 1, -8, -7, -7, -5, 9, 9, 8, 12, -50, -13, -50, -10, -6, -7, -7, -5, -50, -9, -14, -8, -2, -6, -2, -3,
    -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, 33, -50, 30, -50, -50, -50, -50, -50, 29, -50, -50, -50, -50, -50, -50, -50,
    -14, -6, -14, -7, -14, -11, -13, -12, -14, -9, -14, -10, -13, -11, -12, -12, -8, 2, -8, 0, -15, -11, -13, -13, -11, -7, -11, -7, -10, -7, -10, -8, -15, -10, -15, -12, -14, -12, -11, -12, -18, -12, -16, -14, -14, -10, -13, -13, -50, 15, -50, 13, -9, -5, -8, -7, -50, -1, -7, -3, -8, 3, -9, 2,
    -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, 30, -50, 35, -50, -50, -50, -50, -50, 28, -50, -50, -50, -50, -50, -50, -50,
    -13, -7, -14, -5, -12, -11, -11, -10, -13, -9, -12, -7, -11, -11, -11, -10, -7, 0, -8, 2, -14, -12, -13, -11, -9, -8, -11, -5, -9, -8, -9, -7, -13, -11, -14, -8, -12, -12, -12, -10, -15, -13, -14, -11, -12, -11, -13, -10, -50, 13, -50, 15, -8, -7, -8, -4, -50, -3, -7, 0, -7, 2, -8, 3,
    -7, -6, -8, -5, 2, 0, 1, 1, -8, -2, -10, -1, -5, -10, -5, -9, -4, -7, -6, -7, 2, 0, 0, 0, -8, -11, -9, -9, -4, -9, -7, -8, -7, -9, -8, -9, 2, 0, 1, 0, -5, -6, -5, -6, -4, -7, -6, -6, -50, -9, -50, -8, 12, 9, 11, 9, -50, -4, -8, -3, 0, -8, -2, -7,
    -10, -4, -10, -5, 0, 2, 0, 0, -10, 0, -10, -2, -8, -9, -8, -9, -6, -5, -7, -7, -1, 2, 0, 0, -10, -7, -9, -8, -8, -7, -9, -9, -10, -8, -10, -10, 0, 2, 0, 0, -7, -5, -7, -6, -7, -5, -8, -7, -50, -5, -50, -7, 9, 12, 10, 10, -50, -1, -11, -3, -5, -4, -5, -6,
    -8, -6, -8, -5, 1, 0, 2, 1, -9, -1, -8, -1, -7, -10, -4, -9, -5, -7, -5, -7, 0, 0, 1, 0, -9, -8, -7, -8, -5, -9, -6, -8, -8, -9, -8, -9, 1, 0, 3, 0, -6, -5, -4, -5, -5, -7, -5, -7, -50, -8, -50, -8, 11, 10, 13, 10, -50, -4, -6, -3, -3, -8, -1, -6,
    -9, -5, -10, -4, 1, 0, 0, 2, -9, -1, -10, 0, -8, -10, -7, -8, -6, -6, -7, -5, 0, 0, 0, 2, -11, -10, -10, -7, -8, -9, -9, -6, -9, -9, -10, -7, 0, 0, 0, 2, -7, -6, -7, -5, -6, -7, -7, -5, -50, -7, -50, -4, 9, 10, 10, 12, -50, -3, -11, 0, -4, -6, -4, -3,
    -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, -50, 29, -50, 28, -50, -50, -50, -50, -50, 33, -50, -50, -50, -50, -50, -50, -50,
    -13, -7, -13, -8, -9, -6, -8, -7, -9, 0, -9, -2, -12, -10, -12, -11, -9, -4, -10, -5, -13, -10, -11, -10, -6, -1, -6, -3, -12, -8, -12, -9, -16, -12, -15, -14, -8, -6, -7, -7, -9, -4, -8, -6, -11, -8, -11, -9, -50, -1, -50, -3, -4, -1, -4, -3, -50, 16, -5, 14, -10, -4, -10, -5,
    -13, -16, -11, -15, -13, -15, -12, -16, -7, -13, -4, -12, -14, -16, -10, -15, -9, -11, -7, -9, -14, -17, -11, -17, -5, -7, -2, -7, -9, -11, -8, -11, -17, -20, -15, -19, -14, -14, -12, -16, -11, -12, -6, -13, -14, -15, -12, -14, -50, -7, -50, -7, -8, -11, -6, -11, -50, -5, 18, -4, -8, -8, -4, -7,
    -12, -8, -13, -5, -7, -7, -7, -5, -8, -2, -9, 0, -11, -11, -11, -9, -8, -5, -10, -3, -11, -11, -10, -8, -6, -4, -6, -1, -10, -9, -11, -8, -14, -13, -15, -11, -7, -7, -7, -5, -8, -5, -8, -3, -9, -9, -10, -8, -50, -3, -50, 0, -3, -3, -3, 0, -50, 14, -4, 16, -8, -5, -8, -3,
    -10, -13, -12, -11, -4, -8, -6, -7, -10, -11, -11, -11, 2, -1, 0, -1, -6, -9, -7, -8, -5, -8, -6, -8, -10, -11, -9, -10, 9, 6, 7, 6, -12, -16, -13, -16, -5, -8, -7, -7, -12, -14, -12, -13, 1, -3, -1, -2, -50, -8, -50, -7, 0, -5, -3, -4, -50, -10, -8, -8, 13, 0, 9, 0,
    -18, -12, -19, -14, -12, -9, -11, -10, -17, -11, -18, -13, -6, -3, -6, -5, -12, -5, -13, -7, -13, -9, -12, -12, -16, -11, -15, -12, -3, 0, -3, -1, -18, -16, -18, -19, -12, -9, -10, -11, -17, -14, -16, -16, -7, -3, -7, -6, -50, 3, -50, 2, -8, -4, -8, -6, -50, -4, -8, -5, 0, 14, -1, 11,
    -11, -13, -11, -12, -6, -8, -5, -7, -11, -12, -11, -11, 0, -2, 2, -1, -6, -9, -6, -8, -6, -8, -5, -7, -10, -11, -8, -9, 8, 6, 8, 6, -12, -17, -12, -15, -6, -8, -5, -7, -13, -13, -10, -13, 0, -3, 0, -2, -50, -9, -50, -8, -2, -5, -1, -4, -50, -10, -4, -8, 9, -1, 11, 0,
    -17, -13, -17, -11, -11, -11, -10, -9, -16, -12, -16, -11, -5, -5, -6, -3, -11, -7, -13, -5, -12, -11, -11, -9, -15, -15, -14, -11, -2, -1, -2, 0, -17, -18, -18, -15, -10, -11, -10, -9, -15, -14, -15, -13, -6, -5, -7, -3, -50, 2, -50, 3, -7, -6, -6, -3, -50, -5, -7, -3, 0, 11, 0, 14
  ]);

  var CODON_MATRIX = {
    name: "codon",
    size: 64,
    scores: CODON_SCORES
  };

  function safeAdd(score, delta) {
    if (score <= NEG_INF / 2) {
      return NEG_INF;
    }
    return score + delta;
  }

  function safeCombine(left, right) {
    if (left <= NEG_INF / 2 || right <= NEG_INF / 2) {
      return NEG_INF;
    }
    return left + right;
  }

  function encodeCodonSequence(sequence) {
    var upper = sequence.toUpperCase().replace(/U/g, "T");
    var numCodons = upper.length / 3;
    var encoded = new Uint8Array(numCodons);
    var i;
    var codon;
    var idx;

    if (upper.length % 3 !== 0) {
      throw "Sequence length must be divisible by 3.";
    }

    for (i = 0; i < numCodons; i++) {
      codon = upper.substr(i * 3, 3);
      idx = CODON_LOOKUP[codon];
      if (idx === undefined) {
        throw "Unsupported codon '" + codon + "'.";
      }
      encoded[i] = idx;
    }

    return encoded;
  }

  function initializeRow(targetM, targetX, targetY, startMode, gapOpen, gapExtend) {
    var j;

    for (j = 0; j < targetM.length; j++) {
      targetM[j] = NEG_INF;
      targetX[j] = NEG_INF;
      targetY[j] = NEG_INF;
    }

    if (startMode === MODE_START || startMode === STATE_M) {
      targetM[0] = 0;
    } else if (startMode === STATE_X) {
      targetX[0] = 0;
    } else if (startMode === STATE_Y) {
      targetY[0] = 0;
    } else if (startMode === MODE_ANY) {
      targetM[0] = 0;
      targetX[0] = 0;
      targetY[0] = 0;
    } else {
      throw "Unknown alignment state.";
    }

    for (j = 1; j < targetM.length; j++) {
      targetY[j] = Math.max(
        safeAdd(targetM[j - 1], -gapOpen),
        safeAdd(targetY[j - 1], -gapExtend)
      );
    }
  }

  function computeForwardRowScores(seqOne, seqTwo, startMode, gapOpen, gapExtend) {
    var n = seqTwo.length;
    var prevM = new Int32Array(n + 1);
    var prevX = new Int32Array(n + 1);
    var prevY = new Int32Array(n + 1);
    var currM = new Int32Array(n + 1);
    var currX = new Int32Array(n + 1);
    var currY = new Int32Array(n + 1);
    var i;
    var j;
    var rowOffset;
    var residueOne;
    var tmp;

    initializeRow(prevM, prevX, prevY, startMode, gapOpen, gapExtend);

    for (i = 1; i <= seqOne.length; i++) {
      for (j = 0; j <= n; j++) {
        currM[j] = NEG_INF;
        currX[j] = NEG_INF;
        currY[j] = NEG_INF;
      }

      currX[0] = Math.max(
        safeAdd(prevM[0], -gapOpen),
        safeAdd(prevX[0], -gapExtend)
      );

      residueOne = seqOne[i - 1];
      rowOffset = residueOne * 64;

      for (j = 1; j <= n; j++) {
        currX[j] = Math.max(
          safeAdd(prevM[j], -gapOpen),
          safeAdd(prevX[j], -gapExtend)
        );
        currY[j] = Math.max(
          safeAdd(currM[j - 1], -gapOpen),
          safeAdd(currY[j - 1], -gapExtend)
        );

        tmp = prevM[j - 1];
        if (prevX[j - 1] > tmp) { tmp = prevX[j - 1]; }
        if (prevY[j - 1] > tmp) { tmp = prevY[j - 1]; }

        currM[j] = safeAdd(tmp, CODON_SCORES[rowOffset + seqTwo[j - 1]]);
      }

      tmp = prevM; prevM = currM; currM = tmp;
      tmp = prevX; prevX = currX; currX = tmp;
      tmp = prevY; prevY = currY; currY = tmp;
    }

    return { m: prevM, x: prevX, y: prevY };
  }

  function computeReverseRowScores(seqOne, seqTwo, endMode, gapOpen, gapExtend) {
    var n = seqTwo.length;
    var prevM = new Int32Array(n + 1);
    var prevX = new Int32Array(n + 1);
    var prevY = new Int32Array(n + 1);
    var currM = new Int32Array(n + 1);
    var currX = new Int32Array(n + 1);
    var currY = new Int32Array(n + 1);
    var i;
    var j;
    var rowOffset;
    var residueOne;
    var tmp;

    initializeRow(prevM, prevX, prevY, endMode, gapOpen, gapExtend);

    for (i = 1; i <= seqOne.length; i++) {
      for (j = 0; j <= n; j++) {
        currM[j] = NEG_INF;
        currX[j] = NEG_INF;
        currY[j] = NEG_INF;
      }

      currX[0] = Math.max(
        safeAdd(prevM[0], -gapOpen),
        safeAdd(prevX[0], -gapExtend)
      );

      residueOne = seqOne[seqOne.length - i];
      rowOffset = residueOne * 64;

      for (j = 1; j <= n; j++) {
        currX[j] = Math.max(
          safeAdd(prevM[j], -gapOpen),
          safeAdd(prevX[j], -gapExtend)
        );
        currY[j] = Math.max(
          safeAdd(currM[j - 1], -gapOpen),
          safeAdd(currY[j - 1], -gapExtend)
        );

        tmp = prevM[j - 1];
        if (prevX[j - 1] > tmp) { tmp = prevX[j - 1]; }
        if (prevY[j - 1] > tmp) { tmp = prevY[j - 1]; }

        currM[j] = safeAdd(
          tmp,
          CODON_SCORES[rowOffset + seqTwo[seqTwo.length - j]]
        );
      }

      tmp = prevM; prevM = currM; currM = tmp;
      tmp = prevX; prevX = currX; currX = tmp;
      tmp = prevY; prevY = currY; currY = tmp;
    }

    return { m: prevM, x: prevX, y: prevY };
  }

  function chooseBestEndState(scoreM, scoreX, scoreY, endMode) {
    if (endMode === STATE_M) { return { state: STATE_M, score: scoreM }; }
    if (endMode === STATE_X) { return { state: STATE_X, score: scoreX }; }
    if (endMode === STATE_Y) { return { state: STATE_Y, score: scoreY }; }

    if (scoreX > scoreM) {
      return scoreY > scoreX
        ? { state: STATE_Y, score: scoreY }
        : { state: STATE_X, score: scoreX };
    }
    if (scoreY > scoreM) { return { state: STATE_Y, score: scoreY }; }
    return { state: STATE_M, score: scoreM };
  }

  function solveSmall(seqOne, seqTwo, startMode, endMode, gapOpen, gapExtend) {
    var rows = seqOne.length + 1;
    var cols = seqTwo.length + 1;
    var total = rows * cols;
    var scoresM = new Int32Array(total);
    var scoresX = new Int32Array(total);
    var scoresY = new Int32Array(total);
    var traceM = new Uint8Array(total);
    var traceX = new Uint8Array(total);
    var traceY = new Uint8Array(total);
    var i;
    var j;
    var idx;
    var diagIndex;
    var upIndex;
    var leftIndex;
    var scoreFromM;
    var scoreFromX;
    var scoreFromY;
    var best;
    var rowOffset;
    var state;
    var alignedOne = [];
    var alignedTwo = [];

    for (idx = 0; idx < total; idx++) {
      scoresM[idx] = NEG_INF;
      scoresX[idx] = NEG_INF;
      scoresY[idx] = NEG_INF;
      traceM[idx] = 255;
      traceX[idx] = 255;
      traceY[idx] = 255;
    }

    if (startMode === MODE_START || startMode === STATE_M) {
      scoresM[0] = 0;
    } else if (startMode === STATE_X) {
      scoresX[0] = 0;
    } else if (startMode === STATE_Y) {
      scoresY[0] = 0;
    } else if (startMode === MODE_ANY) {
      scoresM[0] = 0;
      scoresX[0] = 0;
      scoresY[0] = 0;
    }

    for (j = 1; j < cols; j++) {
      idx = j;
      leftIndex = j - 1;
      scoreFromM = safeAdd(scoresM[leftIndex], -gapOpen);
      scoreFromY = safeAdd(scoresY[leftIndex], -gapExtend);

      if (scoreFromM >= scoreFromY) {
        scoresY[idx] = scoreFromM;
        traceY[idx] = STATE_M;
      } else {
        scoresY[idx] = scoreFromY;
        traceY[idx] = STATE_Y;
      }
    }

    for (i = 1; i < rows; i++) {
      idx = i * cols;
      upIndex = idx - cols;
      scoreFromM = safeAdd(scoresM[upIndex], -gapOpen);
      scoreFromX = safeAdd(scoresX[upIndex], -gapExtend);

      if (scoreFromM >= scoreFromX) {
        scoresX[idx] = scoreFromM;
        traceX[idx] = STATE_M;
      } else {
        scoresX[idx] = scoreFromX;
        traceX[idx] = STATE_X;
      }

      rowOffset = seqOne[i - 1] * 64;

      for (j = 1; j < cols; j++) {
        idx = i * cols + j;
        diagIndex = idx - cols - 1;
        upIndex = idx - cols;
        leftIndex = idx - 1;

        scoreFromM = safeAdd(scoresM[upIndex], -gapOpen);
        scoreFromX = safeAdd(scoresX[upIndex], -gapExtend);
        if (scoreFromM >= scoreFromX) {
          scoresX[idx] = scoreFromM;
          traceX[idx] = STATE_M;
        } else {
          scoresX[idx] = scoreFromX;
          traceX[idx] = STATE_X;
        }

        scoreFromM = safeAdd(scoresM[leftIndex], -gapOpen);
        scoreFromY = safeAdd(scoresY[leftIndex], -gapExtend);
        if (scoreFromM >= scoreFromY) {
          scoresY[idx] = scoreFromM;
          traceY[idx] = STATE_M;
        } else {
          scoresY[idx] = scoreFromY;
          traceY[idx] = STATE_Y;
        }

        best = scoresM[diagIndex];
        state = STATE_M;
        if (scoresX[diagIndex] > best) { best = scoresX[diagIndex]; state = STATE_X; }
        if (scoresY[diagIndex] > best) { best = scoresY[diagIndex]; state = STATE_Y; }

        scoresM[idx] = safeAdd(best, CODON_SCORES[rowOffset + seqTwo[j - 1]]);
        traceM[idx] = state;
      }
    }

    idx = total - 1;
    best = chooseBestEndState(scoresM[idx], scoresX[idx], scoresY[idx], endMode);
    state = best.state;
    i = seqOne.length;
    j = seqTwo.length;

    while (i > 0 || j > 0) {
      idx = i * cols + j;

      if (state === STATE_M) {
        alignedOne.push(CODON_ORDER[seqOne[i - 1]]);
        alignedTwo.push(CODON_ORDER[seqTwo[j - 1]]);
        state = traceM[idx];
        i--;
        j--;
      } else if (state === STATE_X) {
        alignedOne.push(CODON_ORDER[seqOne[i - 1]]);
        alignedTwo.push("---");
        state = traceX[idx];
        i--;
      } else {
        alignedOne.push("---");
        alignedTwo.push(CODON_ORDER[seqTwo[j - 1]]);
        state = traceY[idx];
        j--;
      }
    }

    alignedOne.reverse();
    alignedTwo.reverse();

    return {
      alignedOne: alignedOne.join(""),
      alignedTwo: alignedTwo.join(""),
      score: best.score
    };
  }

  function solveRecursive(seqOne, seqTwo, startMode, endMode, gapOpen, gapExtend, chunksOne, chunksTwo) {
    var m = seqOne.length;
    var n = seqTwo.length;
    var baseResult;
    var midpoint;
    var forward;
    var reverse;
    var j;
    var scoreM;
    var scoreX;
    var scoreY;
    var bestState = STATE_M;
    var bestSplit = 0;
    var bestScore = NEG_INF;

    if (m === 0 || n === 0 || m * n <= SMALL_PROBLEM_SIZE) {
      baseResult = solveSmall(seqOne, seqTwo, startMode, endMode, gapOpen, gapExtend);
      chunksOne.push(baseResult.alignedOne);
      chunksTwo.push(baseResult.alignedTwo);
      return baseResult.score;
    }

    midpoint = m >>> 1;
    forward = computeForwardRowScores(
      seqOne.subarray(0, midpoint), seqTwo, startMode, gapOpen, gapExtend
    );
    reverse = computeReverseRowScores(
      seqOne.subarray(midpoint), seqTwo, endMode, gapOpen, gapExtend
    );

    for (j = 0; j <= n; j++) {
      scoreM = safeCombine(forward.m[j], reverse.m[n - j]);
      if (scoreM > bestScore) { bestScore = scoreM; bestState = STATE_M; bestSplit = j; }

      scoreX = safeCombine(forward.x[j], reverse.x[n - j]);
      if (scoreX > bestScore) { bestScore = scoreX; bestState = STATE_X; bestSplit = j; }

      scoreY = safeCombine(forward.y[j], reverse.y[n - j]);
      if (scoreY > bestScore) { bestScore = scoreY; bestState = STATE_Y; bestSplit = j; }
    }

    solveRecursive(
      seqOne.subarray(0, midpoint), seqTwo.subarray(0, bestSplit),
      startMode, bestState, gapOpen, gapExtend, chunksOne, chunksTwo
    );
    solveRecursive(
      seqOne.subarray(midpoint), seqTwo.subarray(bestSplit),
      bestState, endMode, gapOpen, gapExtend, chunksOne, chunksTwo
    );

    return bestScore;
  }

  function alignCodon(sequenceOne, sequenceTwo, gapOpen, gapExtend) {
    var encodedOne;
    var encodedTwo;
    var chunksOne = [];
    var chunksTwo = [];
    var score;

    gapOpen = parseInt(gapOpen, 10);
    gapExtend = parseInt(gapExtend, 10);

    if (isNaN(gapOpen) || isNaN(gapExtend) || gapOpen < 0 || gapExtend < 0) {
      throw "Gap opening and extension penalties must be non-negative integers.";
    }

    encodedOne = encodeCodonSequence(sequenceOne);
    encodedTwo = encodeCodonSequence(sequenceTwo);

    score = solveRecursive(
      encodedOne, encodedTwo,
      MODE_START, MODE_ANY,
      gapOpen, gapExtend,
      chunksOne, chunksTwo
    );

    return {
      alignedM: chunksOne.join(""),
      alignedN: chunksTwo.join(""),
      score: score
    };
  }

  global.PairwiseCodonAlignment = {
    align: alignCodon
  };

  if (typeof module !== "undefined" && module.exports) {
    module.exports = global.PairwiseCodonAlignment;
  }
})(this);
