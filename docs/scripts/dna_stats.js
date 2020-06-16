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

function dnaStats(theDocument) {
  var newDna = "";
  var title = "";
  var maxInput = 500000000;

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

  var itemsToCheck = [
    "/g/ (g)1",
    "/a/ (a)1",
    "/t/ (t)1",
    "/c/ (c)1",
    "/n/ (n)1",
    "/u/ (u)1",
    "/r/ (r)1",
    "/y/ (y)1",
    "/s/ (s)1",
    "/w/ (w)1",
    "/k/ (k)1",
    "/m/ (m)1",
    "/b/ (b)1",
    "/d/ (d)1",
    "/h/ (h)1",
    "/v/ (v)1",
    "/g(?=g)/ (gg)2",
    "/g(?=a)/ (ga)2",
    "/g(?=t)/ (gt)2",
    "/g(?=c)/ (gc)2",
    "/g(?=n)/ (gn)2",
    "/a(?=g)/ (ag)2",
    "/a(?=a)/ (aa)2",
    "/a(?=t)/ (at)2",
    "/a(?=c)/ (ac)2",
    "/a(?=n)/ (an)2",
    "/t(?=g)/ (tg)2",
    "/t(?=a)/ (ta)2",
    "/t(?=t)/ (tt)2",
    "/t(?=c)/ (tc)2",
    "/t(?=n)/ (tn)2",
    "/c(?=g)/ (cg)2",
    "/c(?=a)/ (ca)2",
    "/c(?=t)/ (ct)2",
    "/c(?=c)/ (cc)2",
    "/c(?=n)/ (cn)2",
    "/n(?=g)/ (ng)2",
    "/n(?=a)/ (na)2",
    "/n(?=t)/ (nt)2",
    "/n(?=c)/ (nc)2",
    "/n(?=n)/ (nn)2",
    "/g|c/ (g,c)1",
    "/a|t/ (a,t)1",
    "/r|y|s|w|k/ (r,y,s,w,k)1",
    "/b|h|d|v|n/ (b,h,d,v,n)1",
    "/r|y|s|w|k|m|b|d|h|v|n/ (r,y,s,w,k,m,b,d,h,v,n)1",
  ];

  openWindow("DNA Stats");

  var arrayOfFasta = getArrayOfFasta(theDocument.forms[0].elements[0].value);

  for (var i = 0; i < arrayOfFasta.length; i++) {
    newDna = getSequenceFromFasta(arrayOfFasta[i]);
    title = getTitleFromFasta(arrayOfFasta[i]);

    newDna = removeNonDna(newDna);

    outputWindow.document.write(getInfoFromTitleAndSequence(title, newDna));

    writeSequenceStats(newDna, itemsToCheck);

    outputWindow.document.write("<br />\n<br />\n");
  }

  closeWindow();
  return true;
}
