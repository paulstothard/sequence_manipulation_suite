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

function getRestrictionSiteString(type) {
  if (type.toLowerCase() == "standard") {
    return (
      "/aggcct/ (AatI agg|cct)3," +
      "/gacgtc/ (AatII gacgt|c)1," +
      "/tgcgca/ (Acc16I tgc|gca)3," +
      "/cgcg/ (AccII cg|cg)2," +
      "/tccgga/ (AccIII t|ccgga)5," +
      "/aacgtt/ (AclI aa|cgtt)4," +
      "/cacgtg/ (AcvI cac|gtg)3," +
      "/gtac/ (AfaI gt|ac)2," +
      "/agcgct/ (AfeI agc|gct)3," +
      "/cttaag/ (AflII c|ttaag)5," +
      "/accggt/ (AgeI a|ccggt)5," +
      "/actagt/ (AhlI a|ctagt)5," +
      "/gtgcac/ (Alw441 g|tgcac)5," +
      "/agct/ (AluI ag|ct)2," +
      "/agcgct/ (Aor51HI agc|gct)3," +
      "/gggccc/ (ApaI gggcc|c)1," +
      "/gtgcac/ (ApaLI g|tgcac)5," +
      "/ggcgcgcc/ (AscI gg|cgcgcc)6," +
      "/attaat/ (AseI at|taat)4," +
      "/ggtacc/ (Asp718I g|gtacc)5," +
      "/ttcgaa/ (AsuII tt|cgaa)4," +
      "/c[cty]cg[agr]g/ (AvaI c|ycgrg)5," +
      "/tgcgca/ (AviII tgc|gca)3," +
      "/cctagg/ (AvrII c|ctagg)5," +
      "/tggcca/ (BalI tgg|cca)3," +
      "/ggatcc/ (BamHI g|gatcc)5," +
      "/atcgat/ (BanIII at|cgat)4," +
      "/ggcgcc/ (BbeI ggcgc|c)1," +
      "/cacgtg/ (BbrPI cac|gtg)3," +
      "/gcatgc/ (BbuI gcatg|c)1," +
      "/actagt/ (BcuI a|ctagt)5," +
      "/tgatca/ (BclI t|gatca)5," +
      "/ctag/ (BfaI c|tag)3," +
      "/cttaag/ (BfrI c|ttaag)5," +
      "/atgcat/ (BfrBI atg|cat)3," +
      "/agatct/ (BglII a|gatct)5," +
      "/cctagg/ (BlnI c|ctagg)5," +
      "/atcgat/ (BseCI at|cgat)4," +
      "/gcgcgc/ (BsePI g|cgcgc)5," +
      "/cggccg/ (BseX3I c|ggccg)5," +
      "/accggt/ (BshTI a|ccggt)5," +
      "/tgtaca/ (Bsp1407I t|gtaca)5," +
      "/ccatgg/ (Bsp19I c|catgg)5," +
      "/atcgat/ (BspDI at|cgat)4," +
      "/tccgga/ (BspEI t|ccgga)5," +
      "/tgtaca/ (BsrGI t|gtaca)5," +
      "/gcgcgc/ (BssHII g|cgcgc)5," +
      "/cgcg/ (BstUI cg|cg)2," +
      "/atcgat/ (ClaI at|cgat)4," +
      "/gatc/ (DpnII |gatc)4," +
      "/tttaaa/ (DraI ttt|aaa)3," +
      "/cggccg/ (EagI c|ggccg)5," +
      "/gaattc/ (EcoRI g|aattc)5," +
      "/gatatc/ (EcoRV gat|atc)3," +
      "/ggcgcc/ (EgeI ggc|gcc)3," +
      "/ggccggcc/ (FseI ggccgg|cc)2," +
      "/tgcgca/ (FspI tgc|gca)3," +
      "/ggcc/ (HaeIII gg|cc)2," +
      "/gt[cty][agr]ac/ (HincII gty|rac)3," +
      "/aagctt/ (HindIII a|agctt)5," +
      "/ga[acgturyswkmbdhvn]tc/ (HinfI g|antc)4," +
      "/gttaac/ (HpaI gtt|aac)3," +
      "/ccgg/ (HpaII c|cgg)3," +
      "/ggcgcc/ (KasI g|gcgcc)5," +
      "/ggtacc/ (KpnI ggtac|c)1," +
      "/[acgturyswkmbdhvn]gatc[acgturyswkmbdhvn]/ (MboI |gatc)5," +
      "/caattg/ (MfeI c|aattg)5," +
      "/acgcgt/ (MluI a|cgcgt)5," +
      "/tggcca/ (MscI tgg|cca)3," +
      "/ttaa/ (MseI t|taa)3," +
      "/ccgg/ (MspI c|cgg)3," +
      "/gccggc/ (NaeI gcc|ggc)3," +
      "/ggcgcc/ (NarI gg|cgcc)4," +
      "/ccatgg/ (NcoI c|catgg)5," +
      "/catatg/ (NdeI ca|tatg)4," +
      "/gatc/ (NdeII |gatc)4," +
      "/gccggc/ (NgoMIV g|ccggc)5," +
      "/gctagc/ (NheI g|ctagc)5," +
      "/catg/ (NlaIII catg|)0," +
      "/gcggccgc/ (NotI gc|ggccgc)6," +
      "/tcgcga/ (NruI tcg|cga)3," +
      "/atgcat/ (NsiI atgca|t)1," +
      "/ttaattaa/ (PacI ttaat|taa)3," +
      "/acatgt/ (PciI a|catgt)5," +
      "/ggcc/ (PhoI gg|cc)2," +
      "/gtttaaac/ (PmeI gttt|aaac)4," +
      "/cacgtg/ (PmlI cac|gtg)3," +
      "/ttataa/ (PsiI tta|taa)3," +
      "/ctgcag/ (PstI ctgca|g)1," +
      "/cgatcg/ (PvuI cgat|cg)2," +
      "/cagctg/ (PvuII cag|ctg)3," +
      "/gtac/ (RsaI gt|ac)2," +
      "/gagctc/ (SacI gagct|c)1," +
      "/ccgcgg/ (SacII ccgc|gg)2," +
      "/gtcgac/ (SalI g|tcgac)5," +
      "/cctgcagg/ (SbfI cctgca|gg)2," +
      "/agtact/ (ScaI agt|act)3," +
      "/ggcgcc/ (SfoI ggc|gcc)3," +
      "/cccggg/ (SmaI ccc|ggg)3," +
      "/tacgta/ (SnaBI tac|gta)3," +
      "/actagt/ (SpeI a|ctagt)5," +
      "/gcatgc/ (SphI gcatg|c)1," +
      "/aatatt/ (SspI aat|att)3," +
      "/gagctc/ (SstI gagct|c)1," +
      "/ccgcgg/ (SstII ccgc|gg)2," +
      "/aggcct/ (StuI agg|cct)3," +
      "/atttaaat/ (SwaI attt|aaat)4," +
      "/tcga/ (TaqI t|cga)3," +
      "/ctcgag/ (TliI c|tcgag)5," +
      "/attaat/ (VspI at|taat)4," +
      "/tctaga/ (XbaI t|ctaga)5," +
      "/ctcgag/ (XhoI c|tcgag)5," +
      "/cccggg/ (XmaI c|ccggg)5"
    );
  }

  return true;
}
