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

function getGeneticCodeString(type) {
  //  The Standard Code (transl_table=1)
  //    AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  //  Starts = ---M---------------M---------------M----------------------------
  //  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  //  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  //  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

  if (
    type.toLowerCase() == "standard" ||
    type.toLowerCase() == "transl_table=1"
  ) {
    return (
      "/gc[acgturyswkmbdhvn]/=A," +
      "/[tu]g[ctuy]/=C," +
      "/ga[tcuy]/=D," +
      "/ga[agr]/=E," +
      "/[tu][tu][tcuy]/=F," +
      "/gg[acgturyswkmbdhvn]/=G," +
      "/ca[tcuy]/=H," +
      "/a[tu][atcuwmhy]/=I," +
      "/aa[agr]/=K," +
      "/c[tu][acgturyswkmbdhvn]|[tu][tu][agr]|[ctuy][tu][agr]/=L," +
      "/a[tu]g/=M," +
      "/aa[tucy]/=N," +
      "/cc[acgturyswkmbdhvn]/=P," +
      "/ca[agr]/=Q," +
      "/cg[acgturyswkmbdhvn]|ag[agr]|[cam]g[agr]/=R," +
      "/[tu]c[acgturyswkmbdhvn]|ag[ct]/=S," +
      "/ac[acgturyswkmbdhvn]/=T," +
      "/g[tu][acgturyswkmbdhvn]/=V," +
      "/[tu]gg/=W," +
      "/[tu]a[ctuy]/=Y," +
      "/[tu]a[agr]|[tu]ga|[tu][agr]a/=*"
    );
  }

  //  The Vertebrate Mitochondrial Code (transl_table=2)
  //Standard = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  //    AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG
  //  Starts = --------------------------------MMMM---------------M------------
  //  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  //  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  //  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

  if (type.toLowerCase() == "transl_table=2") {
    return (
      "/gc[acgturyswkmbdhvn]/=A," +
      "/[tu]g[ctuy]/=C," +
      "/ga[tcuy]/=D," +
      "/ga[agr]/=E," +
      "/[tu][tu][tcuy]/=F," +
      "/gg[acgturyswkmbdhvn]/=G," +
      "/ca[tcuy]/=H," +
      "/a[tu][tcuy]/=I," +
      "/aa[agr]/=K," +
      "/c[tu][acgturyswkmbdhvn]|[tu][tu][agr]|[ctuy][tu][agr]/=L," +
      "/a[tu][agr]/=M," +
      "/aa[tucy]/=N," +
      "/cc[acgturyswkmbdhvn]/=P," +
      "/ca[agr]/=Q," +
      "/cg[acgturyswkmbdhvn]/=R," +
      "/[tu]c[acgturyswkmbdhvn]|ag[ct]/=S," +
      "/ac[acgturyswkmbdhvn]/=T," +
      "/g[tu][acgturyswkmbdhvn]/=V," +
      "/[tu]g[agr]/=W," +
      "/[tu]a[ctuy]/=Y," +
      "/[tu]a[agr]|ag[agr]/=*"
    );
  }

  //  The Yeast Mitochondrial Code (transl_table=3)
  //Standard = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  //    AAs  = FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  //  Starts = ----------------------------------MM----------------------------
  //  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  //  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  //  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

  if (type.toLowerCase() == "transl_table=3") {
    return (
      "/gc[acgturyswkmbdhvn]/=A," +
      "/[tu]g[ctuy]/=C," +
      "/ga[tcuy]/=D," +
      "/ga[agr]/=E," +
      "/[tu][tu][tcuy]/=F," +
      "/gg[acgturyswkmbdhvn]/=G," +
      "/ca[tcuy]/=H," +
      "/a[tu][tcuy]/=I," +
      "/aa[agr]/=K," +
      "/[tu][tu][agr]/=L," +
      "/a[tu][agr]/=M," +
      "/aa[tucy]/=N," +
      "/cc[acgturyswkmbdhvn]/=P," +
      "/ca[agr]/=Q," +
      "/cg[acgturyswkmbdhvn]|ag[agr]|[cam]g[agr]/=R," +
      "/[tu]c[acgturyswkmbdhvn]|ag[ct]/=S," +
      "/ac[acgturyswkmbdhvn]|c[tu][acgturyswkmbdhvn]/=T," +
      "/g[tu][acgturyswkmbdhvn]/=V," +
      "/[tu]g[agr]/=W," +
      "/[tu]a[ctuy]/=Y," +
      "/[tu]a[agr]/=*"
    );
  }

  //  The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code (transl_table=4)
  //Standard = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  //    AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  //  Starts = --MM---------------M------------MMMM---------------M------------
  //  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  //  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  //  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

  if (type.toLowerCase() == "transl_table=4") {
    return (
      "/gc[acgturyswkmbdhvn]/=A," +
      "/[tu]g[ctuy]/=C," +
      "/ga[tcuy]/=D," +
      "/ga[agr]/=E," +
      "/[tu][tu][tcuy]/=F," +
      "/gg[acgturyswkmbdhvn]/=G," +
      "/ca[tcuy]/=H," +
      "/a[tu][atcuwmhy]/=I," +
      "/aa[agr]/=K," +
      "/c[tu][acgturyswkmbdhvn]|[tu][tu][agr]|[ctuy][tu][agr]/=L," +
      "/a[tu]g/=M," +
      "/aa[tucy]/=N," +
      "/cc[acgturyswkmbdhvn]/=P," +
      "/ca[agr]/=Q," +
      "/cg[acgturyswkmbdhvn]|ag[agr]|[cam]g[agr]/=R," +
      "/[tu]c[acgturyswkmbdhvn]|ag[ct]/=S," +
      "/ac[acgturyswkmbdhvn]/=T," +
      "/g[tu][acgturyswkmbdhvn]/=V," +
      "/[tu]g[agr]/=W," +
      "/[tu]a[ctuy]/=Y," +
      "/[tu]a[agr]/=*"
    );
  }

  //  The Invertebrate Mitochondrial Code (transl_table=5)
  //Standard = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  //    AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG
  //  Starts = ---M----------------------------MMMM---------------M------------
  //  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  //  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  //  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

  if (type.toLowerCase() == "transl_table=5") {
    return (
      "/gc[acgturyswkmbdhvn]/=A," +
      "/[tu]g[ctuy]/=C," +
      "/ga[tcuy]/=D," +
      "/ga[agr]/=E," +
      "/[tu][tu][tcuy]/=F," +
      "/gg[acgturyswkmbdhvn]/=G," +
      "/ca[tcuy]/=H," +
      "/a[tu][tcuy]/=I," +
      "/aa[agr]/=K," +
      "/c[tu][acgturyswkmbdhvn]|[tu][tu][agr]|[ctuy][tu][agr]/=L," +
      "/a[tu][agr]/=M," +
      "/aa[tucy]/=N," +
      "/cc[acgturyswkmbdhvn]/=P," +
      "/ca[agr]/=Q," +
      "/cg[acgturyswkmbdhvn]/=R," +
      "/[tu]c[acgturyswkmbdhvn]|ag[acgturyswkmbdhvn]/=S," +
      "/ac[acgturyswkmbdhvn]/=T," +
      "/g[tu][acgturyswkmbdhvn]/=V," +
      "/[tu]g[agr]/=W," +
      "/[tu]a[ctuy]/=Y," +
      "/[tu]a[agr]/=*"
    );
  }

  //  The Ciliate, Dasycladacean and Hexamita Nuclear Code (transl_table=6)
  //Standard = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  //    AAs  = FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  //  Starts = -----------------------------------M----------------------------
  //  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  //  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  //  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

  if (type.toLowerCase() == "transl_table=6") {
    return (
      "/gc[acgturyswkmbdhvn]/=A," +
      "/[tu]g[ctuy]/=C," +
      "/ga[tcuy]/=D," +
      "/ga[agr]/=E," +
      "/[tu][tu][tcuy]/=F," +
      "/gg[acgturyswkmbdhvn]/=G," +
      "/ca[tcuy]/=H," +
      "/a[tu][atcuwmhy]/=I," +
      "/aa[agr]/=K," +
      "/c[tu][acgturyswkmbdhvn]|[tu][tu][agr]|[ctuy][tu][agr]/=L," +
      "/a[tu]g/=M," +
      "/aa[tucy]/=N," +
      "/cc[acgturyswkmbdhvn]/=P," +
      "/ca[agr]|[tu]a[agr]|[tcuy]a[agr]/=Q," +
      "/cg[acgturyswkmbdhvn]|ag[agr]|[cam]g[agr]/=R," +
      "/[tu]c[acgturyswkmbdhvn]|ag[ct]/=S," +
      "/ac[acgturyswkmbdhvn]/=T," +
      "/g[tu][acgturyswkmbdhvn]/=V," +
      "/[tu]gg/=W," +
      "/[tu]a[ctuy]/=Y," +
      "/[tu]ga/=*"
    );
  }

  //  The Echinoderm and Flatworm Mitochondrial Code (transl_table=9)
  //Standard = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  //    AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG
  //  Starts = -----------------------------------M---------------M------------
  //  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  //  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  //  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

  if (type.toLowerCase() == "transl_table=9") {
    return (
      "/gc[acgturyswkmbdhvn]/=A," +
      "/[tu]g[ctuy]/=C," +
      "/ga[tcuy]/=D," +
      "/ga[agr]/=E," +
      "/[tu][tu][tcuy]/=F," +
      "/gg[acgturyswkmbdhvn]/=G," +
      "/ca[tcuy]/=H," +
      "/a[tu][atcuwmhy]/=I," +
      "/aag/=K," +
      "/c[tu][acgturyswkmbdhvn]|[tu][tu][agr]|[ctuy][tu][agr]/=L," +
      "/a[tu]g/=M," +
      "/aa[atcuwmhy]/=N," +
      "/cc[acgturyswkmbdhvn]/=P," +
      "/ca[agr]/=Q," +
      "/cg[acgturyswkmbdhvn]/=R," +
      "/[tu]c[acgturyswkmbdhvn]|ag[acgturyswkmbdhvn]/=S," +
      "/ac[acgturyswkmbdhvn]/=T," +
      "/g[tu][acgturyswkmbdhvn]/=V," +
      "/[tu]g[agr]/=W," +
      "/[tu]a[ctuy]/=Y," +
      "/[tu]a[agr]/=*"
    );
  }

  //  The Euplotid Nuclear Code (transl_table=10)
  //Standard = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  //    AAs  = FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  //  Starts = -----------------------------------M----------------------------
  //  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  //  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  //  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

  if (type.toLowerCase() == "transl_table=10") {
    return (
      "/gc[acgturyswkmbdhvn]/=A," +
      "/[tu]g[atcuwmhy]/=C," +
      "/ga[tcuy]/=D," +
      "/ga[agr]/=E," +
      "/[tu][tu][tcuy]/=F," +
      "/gg[acgturyswkmbdhvn]/=G," +
      "/ca[tcuy]/=H," +
      "/a[tu][atcuwmhy]/=I," +
      "/aa[agr]/=K," +
      "/c[tu][acgturyswkmbdhvn]|[tu][tu][agr]|[ctuy][tu][agr]/=L," +
      "/a[tu]g/=M," +
      "/aa[tucy]/=N," +
      "/cc[acgturyswkmbdhvn]/=P," +
      "/ca[agr]/=Q," +
      "/cg[acgturyswkmbdhvn]|ag[agr]|[cam]g[agr]/=R," +
      "/[tu]c[acgturyswkmbdhvn]|ag[ct]/=S," +
      "/ac[acgturyswkmbdhvn]/=T," +
      "/g[tu][acgturyswkmbdhvn]/=V," +
      "/[tu]gg/=W," +
      "/[tu]a[ctuy]/=Y," +
      "/[tu]a[agr]/=*"
    );
  }

  //  The Bacterial and Plant Plastid Code (transl_table=11)
  //Standard = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  //    AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  //  Starts = ---M---------------M------------MMMM---------------M------------
  //  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  //  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  //  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

  if (type.toLowerCase() == "transl_table=11") {
    return (
      "/gc[acgturyswkmbdhvn]/=A," +
      "/[tu]g[ctuy]/=C," +
      "/ga[tcuy]/=D," +
      "/ga[agr]/=E," +
      "/[tu][tu][tcuy]/=F," +
      "/gg[acgturyswkmbdhvn]/=G," +
      "/ca[tcuy]/=H," +
      "/a[tu][atcuwmhy]/=I," +
      "/aa[agr]/=K," +
      "/c[tu][acgturyswkmbdhvn]|[tu][tu][agr]|[ctuy][tu][agr]/=L," +
      "/a[tu]g/=M," +
      "/aa[tucy]/=N," +
      "/cc[acgturyswkmbdhvn]/=P," +
      "/ca[agr]/=Q," +
      "/cg[acgturyswkmbdhvn]|ag[agr]|[cam]g[agr]/=R," +
      "/[tu]c[acgturyswkmbdhvn]|ag[ct]/=S," +
      "/ac[acgturyswkmbdhvn]/=T," +
      "/g[tu][acgturyswkmbdhvn]/=V," +
      "/[tu]gg/=W," +
      "/[tu]a[ctuy]/=Y," +
      "/[tu]a[agr]|[tu]ga|[tu][agr]a/=*"
    );
  }

  //  The Alternative Yeast Nuclear Code (transl_table=12)
  //Standard = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  //    AAs  = FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  //  Starts = -------------------M---------------M----------------------------
  //  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  //  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  //  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

  if (type.toLowerCase() == "transl_table=12") {
    return (
      "/gc[acgturyswkmbdhvn]/=A," +
      "/[tu]g[ctuy]/=C," +
      "/ga[tcuy]/=D," +
      "/ga[agr]/=E," +
      "/[tu][tu][tcuy]/=F," +
      "/gg[acgturyswkmbdhvn]/=G," +
      "/ca[tcuy]/=H," +
      "/a[tu][atcuwmhy]/=I," +
      "/aa[agr]/=K," +
      "/c[tu][atcuwmhy]|[tu][tu][agr]|[ctuy][tu]a/=L," +
      "/a[tu]g/=M," +
      "/aa[tucy]/=N," +
      "/cc[acgturyswkmbdhvn]/=P," +
      "/ca[agr]/=Q," +
      "/cg[acgturyswkmbdhvn]|ag[agr]|[cam]g[agr]/=R," +
      "/[tu]c[acgturyswkmbdhvn]|ag[ct]|c[tu]g/=S," +
      "/ac[acgturyswkmbdhvn]/=T," +
      "/g[tu][acgturyswkmbdhvn]/=V," +
      "/[tu]gg/=W," +
      "/[tu]a[ctuy]/=Y," +
      "/[tu]a[agr]|[tu]ga|[tu][agr]a/=*"
    );
  }

  //  The Ascidian Mitochondrial Code (transl_table=13)
  //Standard = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  //    AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG
  //  Starts = ---M------------------------------MM---------------M------------
  //  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  //  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  //  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

  if (type.toLowerCase() == "transl_table=13") {
    return (
      "/gc[acgturyswkmbdhvn]/=A," +
      "/[tu]g[ctuy]/=C," +
      "/ga[tcuy]/=D," +
      "/ga[agr]/=E," +
      "/[tu][tu][tcuy]/=F," +
      "/gg[acgturyswkmbdhvn]|ag[agr]|[agr]g[agr]/=G," +
      "/ca[tcuy]/=H," +
      "/a[tu][tcuy]/=I," +
      "/aa[agr]/=K," +
      "/c[tu][acgturyswkmbdhvn]|[tu][tu][agr]|[ctuy][tu][agr]/=L," +
      "/a[tu][agr]/=M," +
      "/aa[tucy]/=N," +
      "/cc[acgturyswkmbdhvn]/=P," +
      "/ca[agr]/=Q," +
      "/cg[acgturyswkmbdhvn]/=R," +
      "/[tu]c[acgturyswkmbdhvn]|ag[ct]/=S," +
      "/ac[acgturyswkmbdhvn]/=T," +
      "/g[tu][acgturyswkmbdhvn]/=V," +
      "/[tu]g[agr]/=W," +
      "/[tu]a[ctuy]/=Y," +
      "/[tu]a[agr]/=*"
    );
  }

  //  The Alternative Flatworm Mitochondrial Code (transl_table=14)
  //Standard = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  //    AAs  = FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG
  //  Starts = -----------------------------------M----------------------------
  //  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  //  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  //  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

  if (type.toLowerCase() == "transl_table=14") {
    return (
      "/gc[acgturyswkmbdhvn]/=A," +
      "/[tu]g[ctuy]/=C," +
      "/ga[tcuy]/=D," +
      "/ga[agr]/=E," +
      "/[tu][tu][tcuy]/=F," +
      "/gg[acgturyswkmbdhvn]/=G," +
      "/ca[tcuy]/=H," +
      "/a[tu][atcuwmhy]/=I," +
      "/aag/=K," +
      "/c[tu][acgturyswkmbdhvn]|[tu][tu][agr]|[ctuy][tu][agr]/=L," +
      "/a[tu]g/=M," +
      "/aa[atcuwmhy]/=N," +
      "/cc[acgturyswkmbdhvn]/=P," +
      "/ca[agr]/=Q," +
      "/cg[acgturyswkmbdhvn]/=R," +
      "/[tu]c[acgturyswkmbdhvn]|ag[acgturyswkmbdhvn]/=S," +
      "/ac[acgturyswkmbdhvn]/=T," +
      "/g[tu][acgturyswkmbdhvn]/=V," +
      "/[tu]g[agr]/=W," +
      "/[tu]a[atcuwmhy]/=Y," +
      "/[tu]ag/=*"
    );
  }

  //  Blepharisma Nuclear Code (transl_table=15)
  //Standard = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  //    AAs  = FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  //  Starts = -----------------------------------M----------------------------
  //  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  //  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  //  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

  if (type.toLowerCase() == "transl_table=15") {
    return (
      "/gc[acgturyswkmbdhvn]/=A," +
      "/[tu]g[ctuy]/=C," +
      "/ga[tcuy]/=D," +
      "/ga[agr]/=E," +
      "/[tu][tu][tcuy]/=F," +
      "/gg[acgturyswkmbdhvn]/=G," +
      "/ca[tcuy]/=H," +
      "/a[tu][atcuwmhy]/=I," +
      "/aa[agr]/=K," +
      "/c[tu][acgturyswkmbdhvn]|[tu][tu][agr]|[ctuy][tu][agr]/=L," +
      "/a[tu]g/=M," +
      "/aa[tucy]/=N," +
      "/cc[acgturyswkmbdhvn]/=P," +
      "/ca[agr]|[tu]ag|[tcuy]ag/=Q," +
      "/cg[acgturyswkmbdhvn]|ag[agr]|[cam]g[agr]/=R," +
      "/[tu]c[acgturyswkmbdhvn]|ag[ct]/=S," +
      "/ac[acgturyswkmbdhvn]/=T," +
      "/g[tu][acgturyswkmbdhvn]/=V," +
      "/[tu]gg/=W," +
      "/[tu]a[ctuy]/=Y," +
      "/[tu][agr]a/=*"
    );
  }

  //  Chlorophycean Mitochondrial Code (transl_table=16)
  //Standard = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  //    AAs  = FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  //  Starts = -----------------------------------M----------------------------
  //  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  //  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  //  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

  if (type.toLowerCase() == "transl_table=16") {
    return (
      "/gc[acgturyswkmbdhvn]/=A," +
      "/[tu]g[ctuy]/=C," +
      "/ga[tcuy]/=D," +
      "/ga[agr]/=E," +
      "/[tu][tu][tcuy]/=F," +
      "/gg[acgturyswkmbdhvn]/=G," +
      "/ca[tcuy]/=H," +
      "/a[tu][atcuwmhy]/=I," +
      "/aa[agr]/=K," +
      "/c[tu][acgturyswkmbdhvn]|[tu][tu][agr]|[ctuy][tu][agr]|[tu]ag|[tu][atuw]g/=L," +
      "/a[tu]g/=M," +
      "/aa[tucy]/=N," +
      "/cc[acgturyswkmbdhvn]/=P," +
      "/ca[agr]/=Q," +
      "/cg[acgturyswkmbdhvn]|ag[agr]|[cam]g[agr]/=R," +
      "/[tu]c[acgturyswkmbdhvn]|ag[ct]/=S," +
      "/ac[acgturyswkmbdhvn]/=T," +
      "/g[tu][acgturyswkmbdhvn]/=V," +
      "/[tu]gg/=W," +
      "/[tu]a[ctuy]/=Y," +
      "/[tu][agr]a/=*"
    );
  }

  //  Trematode Mitochondrial Code (transl_table=21)
  //Standard = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  //    AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG
  //  Starts = -----------------------------------M---------------M------------
  //  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  //  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  //  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

  if (type.toLowerCase() == "transl_table=21") {
    return (
      "/gc[acgturyswkmbdhvn]/=A," +
      "/[tu]g[ctuy]/=C," +
      "/ga[tcuy]/=D," +
      "/ga[agr]/=E," +
      "/[tu][tu][tcuy]/=F," +
      "/gg[acgturyswkmbdhvn]/=G," +
      "/ca[tcuy]/=H," +
      "/a[tu][tcuy]/=I," +
      "/aag/=K," +
      "/c[tu][acgturyswkmbdhvn]|[tu][tu][agr]|[ctuy][tu][agr]/=L," +
      "/a[tu][agr]/=M," +
      "/aa[atcuwmhy]/=N," +
      "/cc[acgturyswkmbdhvn]/=P," +
      "/ca[agr]/=Q," +
      "/cg[acgturyswkmbdhvn]/=R," +
      "/[tu]c[acgturyswkmbdhvn]|ag[acgturyswkmbdhvn]/=S," +
      "/ac[acgturyswkmbdhvn]/=T," +
      "/g[tu][acgturyswkmbdhvn]/=V," +
      "/[tu]g[agr]/=W," +
      "/[tu]a[ctuy]/=Y," +
      "/[tu]a[agr]/=*"
    );
  }

  //  Scenedesmus obliquus mitochondrial Code (transl_table=22)
  //Standard = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  //    AAs  = FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  //  Starts = -----------------------------------M----------------------------
  //  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  //  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  //  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

  if (type.toLowerCase() == "transl_table=22") {
    return (
      "/gc[acgturyswkmbdhvn]/=A," +
      "/[tu]g[ctuy]/=C," +
      "/ga[tcuy]/=D," +
      "/ga[agr]/=E," +
      "/[tu][tu][tcuy]/=F," +
      "/gg[acgturyswkmbdhvn]/=G," +
      "/ca[tcuy]/=H," +
      "/a[tu][atcuwmhy]/=I," +
      "/aa[agr]/=K," +
      "/c[tu][acgturyswkmbdhvn]|[tu][tu][agr]|[ctuy][tu][agr]|[tu]ag|[tu][atuw]g/=L," +
      "/a[tu]g/=M," +
      "/aa[tucy]/=N," +
      "/cc[acgturyswkmbdhvn]/=P," +
      "/ca[agr]/=Q," +
      "/cg[acgturyswkmbdhvn]|ag[agr]|[cam]g[agr]/=R," +
      "/[tu]c[cgtyskb]|ag[ct]/=S," +
      "/ac[acgturyswkmbdhvn]/=T," +
      "/g[tu][acgturyswkmbdhvn]/=V," +
      "/[tu]gg/=W," +
      "/[tu]a[ctuy]/=Y," +
      "/[tu][agcrsmv]a/=*"
    );
  }

  //  Thraustochytrium Mitochondrial Code (transl_table=23)
  //Standard = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  //    AAs  = FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  //  Starts = --------------------------------M--M---------------M------------
  //  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  //  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  //  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

  if (type.toLowerCase() == "transl_table=23") {
    return (
      "/gc[acgturyswkmbdhvn]/=A," +
      "/[tu]g[ctuy]/=C," +
      "/ga[tcuy]/=D," +
      "/ga[agr]/=E," +
      "/[tu][tu][tcuy]/=F," +
      "/gg[acgturyswkmbdhvn]/=G," +
      "/ca[tcuy]/=H," +
      "/a[tu][atcuwmhy]/=I," +
      "/aa[agr]/=K," +
      "/c[tu][acgturyswkmbdhvn]|[ctuy][tu]g/=L," +
      "/a[tu]g/=M," +
      "/aa[tucy]/=N," +
      "/cc[acgturyswkmbdhvn]/=P," +
      "/ca[agr]/=Q," +
      "/cg[acgturyswkmbdhvn]|ag[agr]|[cam]g[agr]/=R," +
      "/[tu]c[acgturyswkmbdhvn]|ag[ct]/=S," +
      "/ac[acgturyswkmbdhvn]/=T," +
      "/g[tu][acgturyswkmbdhvn]/=V," +
      "/[tu]gg/=W," +
      "/[tu]a[ctuy]/=Y," +
      "/[tu]a[agr]|[tu]ga|[tu][agtrwkd]a/=*"
    );
  }

  return true;
}
