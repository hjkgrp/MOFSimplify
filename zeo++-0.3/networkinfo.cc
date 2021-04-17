//#include "network.h"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <set>
#include <algorithm>

#include "networkinfo.h"

using namespace std;

/* This file is used to store information about the radii and weights of
 * atoms. The information can be specified in input files or default values can be used
 * for more common types.
 */

/*
 * Parameters taken from CCDC tables
 * */
static std::map <std::string,double> radTable;
static std::map <std::string,double> covRadTable;
static std::map <std::string,double> massTable;
static std::map <std::string,int> atomicNumberTable;
static std::map <std::string,bool> atomicCharacterTable;
static std::set<string> periodicTable;

bool stripAtomNameInternalFlag = false;

/** Initialize atom name stripping */
void initializeStripAtomNameInternalFlag(bool value){
stripAtomNameInternalFlag = value;
}

/** Fills the radius table with several default values. */
void initializeRadTable(){
//radTable.insert(pair <string,double> ("Symbol",  vdW_Radius));
radTable.insert(pair <string,double> ("H",  1.09));
radTable.insert(pair <string,double> ("D",  1.09));
radTable.insert(pair <string,double> ("He",  1.4));
radTable.insert(pair <string,double> ("Li",  1.82));
radTable.insert(pair <string,double> ("Be",  2));
radTable.insert(pair <string,double> ("B",  2));
radTable.insert(pair <string,double> ("C",  1.7));
radTable.insert(pair <string,double> ("N",  1.55));
radTable.insert(pair <string,double> ("O",  1.52));
radTable.insert(pair <string,double> ("F",  1.47));
radTable.insert(pair <string,double> ("Ne",  1.54));
radTable.insert(pair <string,double> ("Na",  2.27));
radTable.insert(pair <string,double> ("Mg",  1.73));
radTable.insert(pair <string,double> ("Al",  2));
radTable.insert(pair <string,double> ("Si",  2.1));
radTable.insert(pair <string,double> ("P",  1.8));
radTable.insert(pair <string,double> ("S",  1.8));
radTable.insert(pair <string,double> ("Cl",  1.75));
radTable.insert(pair <string,double> ("Ar",  1.88));
radTable.insert(pair <string,double> ("K",  2.75));
radTable.insert(pair <string,double> ("Ca",  2));
radTable.insert(pair <string,double> ("Sc",  2));
radTable.insert(pair <string,double> ("Ti",  2));
radTable.insert(pair <string,double> ("V",  2));
radTable.insert(pair <string,double> ("Cr",  2));
radTable.insert(pair <string,double> ("Mn",  2));
radTable.insert(pair <string,double> ("Fe",  2));
radTable.insert(pair <string,double> ("Co",  2));
radTable.insert(pair <string,double> ("Ni",  1.63));
radTable.insert(pair <string,double> ("Cu",  1.4));
radTable.insert(pair <string,double> ("Zn",  1.39));
radTable.insert(pair <string,double> ("Ga",  1.87));
radTable.insert(pair <string,double> ("Ge",  2));
radTable.insert(pair <string,double> ("As",  1.85));
radTable.insert(pair <string,double> ("Se",  1.9));
radTable.insert(pair <string,double> ("Br",  1.85));
radTable.insert(pair <string,double> ("Kr",  2.02));
radTable.insert(pair <string,double> ("Rb",  2));
radTable.insert(pair <string,double> ("Sr",  2));
radTable.insert(pair <string,double> ("Y",  2));
radTable.insert(pair <string,double> ("Zr",  2));
radTable.insert(pair <string,double> ("Nb",  2));
radTable.insert(pair <string,double> ("Mo",  2));
radTable.insert(pair <string,double> ("Tc",  2));
radTable.insert(pair <string,double> ("Ru",  2));
radTable.insert(pair <string,double> ("Rh",  2));
radTable.insert(pair <string,double> ("Pd",  1.63));
radTable.insert(pair <string,double> ("Ag",  1.72));
radTable.insert(pair <string,double> ("Cd",  1.58));
radTable.insert(pair <string,double> ("In",  1.93));
radTable.insert(pair <string,double> ("Sn",  2.17));
radTable.insert(pair <string,double> ("Sb",  2));
radTable.insert(pair <string,double> ("Te",  2.06));
radTable.insert(pair <string,double> ("I",  1.98));
radTable.insert(pair <string,double> ("Xe",  2.16));
radTable.insert(pair <string,double> ("Cs",  2));
radTable.insert(pair <string,double> ("Ba",  2));
radTable.insert(pair <string,double> ("La",  2));
radTable.insert(pair <string,double> ("Ce",  2));
radTable.insert(pair <string,double> ("Pr",  2));
radTable.insert(pair <string,double> ("Nd",  2));
radTable.insert(pair <string,double> ("Pm",  2));
radTable.insert(pair <string,double> ("Sm",  2));
radTable.insert(pair <string,double> ("Eu",  2));
radTable.insert(pair <string,double> ("Gd",  2));
radTable.insert(pair <string,double> ("Tb",  2));
radTable.insert(pair <string,double> ("Dy",  2));
radTable.insert(pair <string,double> ("Ho",  2));
radTable.insert(pair <string,double> ("Er",  2));
radTable.insert(pair <string,double> ("Tm",  2));
radTable.insert(pair <string,double> ("Yb",  2));
radTable.insert(pair <string,double> ("Lu",  2));
radTable.insert(pair <string,double> ("Hf",  2));
radTable.insert(pair <string,double> ("Ta",  2));
radTable.insert(pair <string,double> ("W",  2));
radTable.insert(pair <string,double> ("Re",  2));
radTable.insert(pair <string,double> ("Os",  2));
radTable.insert(pair <string,double> ("Ir",  2));
radTable.insert(pair <string,double> ("Pt",  1.72));
radTable.insert(pair <string,double> ("Au",  1.66));
radTable.insert(pair <string,double> ("Hg",  1.55));
radTable.insert(pair <string,double> ("Tl",  1.96));
radTable.insert(pair <string,double> ("Pb",  2.02));
radTable.insert(pair <string,double> ("Bi",  2));
radTable.insert(pair <string,double> ("Po",  2));
radTable.insert(pair <string,double> ("At",  2));
radTable.insert(pair <string,double> ("Rn",  2));
radTable.insert(pair <string,double> ("Fr",  2));
radTable.insert(pair <string,double> ("Ra",  2));
radTable.insert(pair <string,double> ("Ac",  2));
radTable.insert(pair <string,double> ("Th",  2));
radTable.insert(pair <string,double> ("Pa",  2));
radTable.insert(pair <string,double> ("U",  1.86));
radTable.insert(pair <string,double> ("Np",  2));
radTable.insert(pair <string,double> ("Pu",  2));
radTable.insert(pair <string,double> ("Am",  2));
radTable.insert(pair <string,double> ("Cm",  2));
radTable.insert(pair <string,double> ("Bk",  2));
radTable.insert(pair <string,double> ("Cf",  2));
radTable.insert(pair <string,double> ("Es",  2));
radTable.insert(pair <string,double> ("Fm",  2));
radTable.insert(pair <string,double> ("Md",  2));
radTable.insert(pair <string,double> ("No",  2));
radTable.insert(pair <string,double> ("Lr",  2));
radTable.insert(pair <string,double> ("Rf",  2));
radTable.insert(pair <string,double> ("Db",  2));
radTable.insert(pair <string,double> ("Sg",  2));
radTable.insert(pair <string,double> ("Bh",  2));
radTable.insert(pair <string,double> ("Hs",  2));
radTable.insert(pair <string,double> ("Mt",  2));
radTable.insert(pair <string,double> ("Ds",  2));
}


/** Fills the covalent radius table with several default values. */
/* 
 * When determining bonding using these radii CCDC advises to use
 * threshold +0.4A
 * */
void initializeCovRadTable(){
//covRadTable.insert(pair <string,double> ("Symbol",  Covalent_Radius));
covRadTable.insert(pair <string,double> ("H",  0.23));
covRadTable.insert(pair <string,double> ("D",  0.23));
covRadTable.insert(pair <string,double> ("He",  1.5));
covRadTable.insert(pair <string,double> ("Li",  1.28));
covRadTable.insert(pair <string,double> ("Be",  0.96));
covRadTable.insert(pair <string,double> ("B",  0.83));
covRadTable.insert(pair <string,double> ("C",  0.68));
covRadTable.insert(pair <string,double> ("N",  0.68));
covRadTable.insert(pair <string,double> ("O",  0.68));
covRadTable.insert(pair <string,double> ("F",  0.64));
covRadTable.insert(pair <string,double> ("Ne",  1.5));
covRadTable.insert(pair <string,double> ("Na",  1.66));
covRadTable.insert(pair <string,double> ("Mg",  1.41));
covRadTable.insert(pair <string,double> ("Al",  1.21));
covRadTable.insert(pair <string,double> ("Si",  1.2));
covRadTable.insert(pair <string,double> ("P",  1.05));
covRadTable.insert(pair <string,double> ("S",  1.02));
covRadTable.insert(pair <string,double> ("Cl",  0.99));
covRadTable.insert(pair <string,double> ("Ar",  1.51));
covRadTable.insert(pair <string,double> ("K",  2.03));
covRadTable.insert(pair <string,double> ("Ca",  1.76));
covRadTable.insert(pair <string,double> ("Sc",  1.7));
covRadTable.insert(pair <string,double> ("Ti",  1.6));
covRadTable.insert(pair <string,double> ("V",  1.53));
covRadTable.insert(pair <string,double> ("Cr",  1.39));
covRadTable.insert(pair <string,double> ("Mn",  1.61));
covRadTable.insert(pair <string,double> ("Fe",  1.52));
covRadTable.insert(pair <string,double> ("Co",  1.26));
covRadTable.insert(pair <string,double> ("Ni",  1.24));
covRadTable.insert(pair <string,double> ("Cu",  1.32));
covRadTable.insert(pair <string,double> ("Zn",  1.22));
covRadTable.insert(pair <string,double> ("Ga",  1.22));
covRadTable.insert(pair <string,double> ("Ge",  1.17));
covRadTable.insert(pair <string,double> ("As",  1.21));
covRadTable.insert(pair <string,double> ("Se",  1.22));
covRadTable.insert(pair <string,double> ("Br",  1.21));
covRadTable.insert(pair <string,double> ("Kr",  1.5));
covRadTable.insert(pair <string,double> ("Rb",  2.2));
covRadTable.insert(pair <string,double> ("Sr",  1.95));
covRadTable.insert(pair <string,double> ("Y",  1.9));
covRadTable.insert(pair <string,double> ("Zr",  1.75));
covRadTable.insert(pair <string,double> ("Nb",  1.64));
covRadTable.insert(pair <string,double> ("Mo",  1.54));
covRadTable.insert(pair <string,double> ("Tc",  1.47));
covRadTable.insert(pair <string,double> ("Ru",  1.46));
covRadTable.insert(pair <string,double> ("Rh",  1.42));
covRadTable.insert(pair <string,double> ("Pd",  1.39));
covRadTable.insert(pair <string,double> ("Ag",  1.45));
covRadTable.insert(pair <string,double> ("Cd",  1.54));
covRadTable.insert(pair <string,double> ("In",  1.42));
covRadTable.insert(pair <string,double> ("Sn",  1.39));
covRadTable.insert(pair <string,double> ("Sb",  1.39));
covRadTable.insert(pair <string,double> ("Te",  1.47));
covRadTable.insert(pair <string,double> ("I",  1.4));
covRadTable.insert(pair <string,double> ("Xe",  1.5));
covRadTable.insert(pair <string,double> ("Cs",  2.44));
covRadTable.insert(pair <string,double> ("Ba",  2.15));
covRadTable.insert(pair <string,double> ("La",  2.07));
covRadTable.insert(pair <string,double> ("Ce",  2.04));
covRadTable.insert(pair <string,double> ("Pr",  2.03));
covRadTable.insert(pair <string,double> ("Nd",  2.01));
covRadTable.insert(pair <string,double> ("Pm",  1.99));
covRadTable.insert(pair <string,double> ("Sm",  1.98));
covRadTable.insert(pair <string,double> ("Eu",  1.98));
covRadTable.insert(pair <string,double> ("Gd",  1.96));
covRadTable.insert(pair <string,double> ("Tb",  1.94));
covRadTable.insert(pair <string,double> ("Dy",  1.92));
covRadTable.insert(pair <string,double> ("Ho",  1.92));
covRadTable.insert(pair <string,double> ("Er",  1.89));
covRadTable.insert(pair <string,double> ("Tm",  1.9));
covRadTable.insert(pair <string,double> ("Yb",  1.87));
covRadTable.insert(pair <string,double> ("Lu",  1.87));
covRadTable.insert(pair <string,double> ("Hf",  1.75));
covRadTable.insert(pair <string,double> ("Ta",  1.7));
covRadTable.insert(pair <string,double> ("W",  1.62));
covRadTable.insert(pair <string,double> ("Re",  1.51));
covRadTable.insert(pair <string,double> ("Os",  1.44));
covRadTable.insert(pair <string,double> ("Ir",  1.41));
covRadTable.insert(pair <string,double> ("Pt",  1.36));
covRadTable.insert(pair <string,double> ("Au",  1.36));
covRadTable.insert(pair <string,double> ("Hg",  1.32));
covRadTable.insert(pair <string,double> ("Tl",  1.45));
covRadTable.insert(pair <string,double> ("Pb",  1.46));
covRadTable.insert(pair <string,double> ("Bi",  1.48));
covRadTable.insert(pair <string,double> ("Po",  1.4));
covRadTable.insert(pair <string,double> ("At",  1.21));
covRadTable.insert(pair <string,double> ("Rn",  1.5));
covRadTable.insert(pair <string,double> ("Fr",  2.6));
covRadTable.insert(pair <string,double> ("Ra",  2.21));
covRadTable.insert(pair <string,double> ("Ac",  2.15));
covRadTable.insert(pair <string,double> ("Th",  2.06));
covRadTable.insert(pair <string,double> ("Pa",  2));
covRadTable.insert(pair <string,double> ("U",  1.96));
covRadTable.insert(pair <string,double> ("Np",  1.9));
covRadTable.insert(pair <string,double> ("Pu",  1.87));
covRadTable.insert(pair <string,double> ("Am",  1.8));
covRadTable.insert(pair <string,double> ("Cm",  1.69));
covRadTable.insert(pair <string,double> ("Bk",  1.54));
covRadTable.insert(pair <string,double> ("Cf",  1.83));
covRadTable.insert(pair <string,double> ("Es",  1.5));
covRadTable.insert(pair <string,double> ("Fm",  1.5));
covRadTable.insert(pair <string,double> ("Md",  1.5));
covRadTable.insert(pair <string,double> ("No",  1.5));
covRadTable.insert(pair <string,double> ("Lr",  1.5));
covRadTable.insert(pair <string,double> ("Rf",  1.5));
covRadTable.insert(pair <string,double> ("Db",  1.5));
covRadTable.insert(pair <string,double> ("Sg",  1.5));
covRadTable.insert(pair <string,double> ("Bh",  1.5));
covRadTable.insert(pair <string,double> ("Hs",  1.5));
covRadTable.insert(pair <string,double> ("Mt",  1.5));
covRadTable.insert(pair <string,double> ("Ds",  1.5));
}

/** Fills the mass table with several default values
 ** in units of g/mole. */
void initializeMassTable(){
//massTable.insert(pair <string,double> ("Symbol",  Atomic_Weight));
massTable.insert(pair <string,double> ("H",  1.008));
massTable.insert(pair <string,double> ("D",  2.0));
massTable.insert(pair <string,double> ("He",  4.003));
massTable.insert(pair <string,double> ("Li",  6.941));
massTable.insert(pair <string,double> ("Be",  9.012));
massTable.insert(pair <string,double> ("B",  10.811));
massTable.insert(pair <string,double> ("C",  12.011));
massTable.insert(pair <string,double> ("N",  14.007));
massTable.insert(pair <string,double> ("O",  15.999));
massTable.insert(pair <string,double> ("F",  18.998));
massTable.insert(pair <string,double> ("Ne",  20.18));
massTable.insert(pair <string,double> ("Na",  22.991));
massTable.insert(pair <string,double> ("Mg",  24.305));
massTable.insert(pair <string,double> ("Al",  26.982));
massTable.insert(pair <string,double> ("Si",  28.086));
massTable.insert(pair <string,double> ("P",  30.974));
massTable.insert(pair <string,double> ("S",  32.066));
massTable.insert(pair <string,double> ("Cl",  35.453));
massTable.insert(pair <string,double> ("Ar",  39.948));
massTable.insert(pair <string,double> ("K",  39.098));
massTable.insert(pair <string,double> ("Ca",  40.078));
massTable.insert(pair <string,double> ("Sc",  44.956));
massTable.insert(pair <string,double> ("Ti",  47.867));
massTable.insert(pair <string,double> ("V",  50.942));
massTable.insert(pair <string,double> ("Cr",  51.996));
massTable.insert(pair <string,double> ("Mn",  54.938));
massTable.insert(pair <string,double> ("Fe",  55.845));
massTable.insert(pair <string,double> ("Co",  58.933));
massTable.insert(pair <string,double> ("Ni",  58.693));
massTable.insert(pair <string,double> ("Cu",  63.546));
massTable.insert(pair <string,double> ("Zn",  65.39));
massTable.insert(pair <string,double> ("Ga",  69.723));
massTable.insert(pair <string,double> ("Ge",  72.61));
massTable.insert(pair <string,double> ("As",  74.922));
massTable.insert(pair <string,double> ("Se",  78.96));
massTable.insert(pair <string,double> ("Br",  79.904));
massTable.insert(pair <string,double> ("Kr",  83.8));
massTable.insert(pair <string,double> ("Rb",  85.468));
massTable.insert(pair <string,double> ("Sr",  87.62));
massTable.insert(pair <string,double> ("Y",  88.906));
massTable.insert(pair <string,double> ("Zr",  91.224));
massTable.insert(pair <string,double> ("Nb",  92.906));
massTable.insert(pair <string,double> ("Mo",  95.94));
massTable.insert(pair <string,double> ("Tc",  98));
massTable.insert(pair <string,double> ("Ru",  101.07));
massTable.insert(pair <string,double> ("Rh",  102.906));
massTable.insert(pair <string,double> ("Pd",  106.42));
massTable.insert(pair <string,double> ("Ag",  107.868));
massTable.insert(pair <string,double> ("Cd",  112.411));
massTable.insert(pair <string,double> ("In",  114.818));
massTable.insert(pair <string,double> ("Sn",  118.71));
massTable.insert(pair <string,double> ("Sb",  121.76));
massTable.insert(pair <string,double> ("Te",  127.6));
massTable.insert(pair <string,double> ("I",  126.904));
massTable.insert(pair <string,double> ("Xe",  131.29));
massTable.insert(pair <string,double> ("Cs",  132.905));
massTable.insert(pair <string,double> ("Ba",  137.327));
massTable.insert(pair <string,double> ("La",  138.906));
massTable.insert(pair <string,double> ("Ce",  140.116));
massTable.insert(pair <string,double> ("Pr",  140.908));
massTable.insert(pair <string,double> ("Nd",  144.24));
massTable.insert(pair <string,double> ("Pm",  145));
massTable.insert(pair <string,double> ("Sm",  150.36));
massTable.insert(pair <string,double> ("Eu",  151.964));
massTable.insert(pair <string,double> ("Gd",  157.25));
massTable.insert(pair <string,double> ("Tb",  158.925));
massTable.insert(pair <string,double> ("Dy",  162.5));
massTable.insert(pair <string,double> ("Ho",  164.93));
massTable.insert(pair <string,double> ("Er",  167.26));
massTable.insert(pair <string,double> ("Tm",  168.934));
massTable.insert(pair <string,double> ("Yb",  173.04));
massTable.insert(pair <string,double> ("Lu",  174.967));
massTable.insert(pair <string,double> ("Hf",  178.49));
massTable.insert(pair <string,double> ("Ta",  180.948));
massTable.insert(pair <string,double> ("W",  183.84));
massTable.insert(pair <string,double> ("Re",  186.207));
massTable.insert(pair <string,double> ("Os",  190.23));
massTable.insert(pair <string,double> ("Ir",  192.217));
massTable.insert(pair <string,double> ("Pt",  195.078));
massTable.insert(pair <string,double> ("Au",  196.967));
massTable.insert(pair <string,double> ("Hg",  200.59));
massTable.insert(pair <string,double> ("Tl",  204.383));
massTable.insert(pair <string,double> ("Pb",  207.2));
massTable.insert(pair <string,double> ("Bi",  208.98));
massTable.insert(pair <string,double> ("Po",  210));
massTable.insert(pair <string,double> ("At",  210));
massTable.insert(pair <string,double> ("Rn",  222));
massTable.insert(pair <string,double> ("Fr",  223));
massTable.insert(pair <string,double> ("Ra",  226));
massTable.insert(pair <string,double> ("Ac",  227));
massTable.insert(pair <string,double> ("Th",  232.038));
massTable.insert(pair <string,double> ("Pa",  231.036));
massTable.insert(pair <string,double> ("U",  238.029));
massTable.insert(pair <string,double> ("Np",  237));
massTable.insert(pair <string,double> ("Pu",  244));
massTable.insert(pair <string,double> ("Am",  243));
massTable.insert(pair <string,double> ("Cm",  247));
massTable.insert(pair <string,double> ("Bk",  247));
massTable.insert(pair <string,double> ("Cf",  251));
massTable.insert(pair <string,double> ("Es",  252));
massTable.insert(pair <string,double> ("Fm",  257));
massTable.insert(pair <string,double> ("Md",  258));
massTable.insert(pair <string,double> ("No",  259));
massTable.insert(pair <string,double> ("Lr",  262));
massTable.insert(pair <string,double> ("Rf",  261));
massTable.insert(pair <string,double> ("Db",  262));
massTable.insert(pair <string,double> ("Sg",  266));
massTable.insert(pair <string,double> ("Bh",  264));
massTable.insert(pair <string,double> ("Hs",  269));
massTable.insert(pair <string,double> ("Mt",  268));
massTable.insert(pair <string,double> ("Ds",  271));
}


/** Fills atom metal/nonmetal information **/
void  initializeAtomCharacterTable(){
  // value set to true if atom is metal

//  atomicCharacterTable.insert(pair <string,bool> ("Symbol",true));
  atomicCharacterTable.insert(pair <string,bool> ("H",false));
  atomicCharacterTable.insert(pair <string,bool> ("D",false));
  atomicCharacterTable.insert(pair <string,bool> ("He",false));
  atomicCharacterTable.insert(pair <string,bool> ("Li",true));
  atomicCharacterTable.insert(pair <string,bool> ("Be",true));
  atomicCharacterTable.insert(pair <string,bool> ("B",false));
  atomicCharacterTable.insert(pair <string,bool> ("C",false));
  atomicCharacterTable.insert(pair <string,bool> ("N",false));
  atomicCharacterTable.insert(pair <string,bool> ("O",false));
  atomicCharacterTable.insert(pair <string,bool> ("F",false));
  atomicCharacterTable.insert(pair <string,bool> ("Ne",false));
  atomicCharacterTable.insert(pair <string,bool> ("Na",true));
  atomicCharacterTable.insert(pair <string,bool> ("Mg",true));
  atomicCharacterTable.insert(pair <string,bool> ("Al",true));
  atomicCharacterTable.insert(pair <string,bool> ("Si",false));
  atomicCharacterTable.insert(pair <string,bool> ("P",false));
  atomicCharacterTable.insert(pair <string,bool> ("S",false));
  atomicCharacterTable.insert(pair <string,bool> ("Cl",false));
  atomicCharacterTable.insert(pair <string,bool> ("Ar",false));
  atomicCharacterTable.insert(pair <string,bool> ("K",true));
  atomicCharacterTable.insert(pair <string,bool> ("Ca",true));
  atomicCharacterTable.insert(pair <string,bool> ("Sc",true));
  atomicCharacterTable.insert(pair <string,bool> ("Ti",true));
  atomicCharacterTable.insert(pair <string,bool> ("V",true));
  atomicCharacterTable.insert(pair <string,bool> ("Cr",true));
  atomicCharacterTable.insert(pair <string,bool> ("Mn",true));
  atomicCharacterTable.insert(pair <string,bool> ("Fe",true));
  atomicCharacterTable.insert(pair <string,bool> ("Co",true));
  atomicCharacterTable.insert(pair <string,bool> ("Ni",true));
  atomicCharacterTable.insert(pair <string,bool> ("Cu",true));
  atomicCharacterTable.insert(pair <string,bool> ("Zn",true));
  atomicCharacterTable.insert(pair <string,bool> ("Ga",true));
  atomicCharacterTable.insert(pair <string,bool> ("Ge",false));
  atomicCharacterTable.insert(pair <string,bool> ("As",false));
  atomicCharacterTable.insert(pair <string,bool> ("Se",false));
  atomicCharacterTable.insert(pair <string,bool> ("Br",false));
  atomicCharacterTable.insert(pair <string,bool> ("Kr",false));
  atomicCharacterTable.insert(pair <string,bool> ("Rb",true));
  atomicCharacterTable.insert(pair <string,bool> ("Sr",true));
  atomicCharacterTable.insert(pair <string,bool> ("Y",true));
  atomicCharacterTable.insert(pair <string,bool> ("Zr",true));
  atomicCharacterTable.insert(pair <string,bool> ("Nb",true));
  atomicCharacterTable.insert(pair <string,bool> ("Mo",true));
  atomicCharacterTable.insert(pair <string,bool> ("Tc",true));
  atomicCharacterTable.insert(pair <string,bool> ("Ru",true));
  atomicCharacterTable.insert(pair <string,bool> ("Rh",true));
  atomicCharacterTable.insert(pair <string,bool> ("Pd",true));
  atomicCharacterTable.insert(pair <string,bool> ("Ag",true));
  atomicCharacterTable.insert(pair <string,bool> ("Cd",true));
  atomicCharacterTable.insert(pair <string,bool> ("In",true));
  atomicCharacterTable.insert(pair <string,bool> ("Sn",true));
  atomicCharacterTable.insert(pair <string,bool> ("Sb",true));
  atomicCharacterTable.insert(pair <string,bool> ("Te",true));
  atomicCharacterTable.insert(pair <string,bool> ("I",false));
  atomicCharacterTable.insert(pair <string,bool> ("Xe",false));
  atomicCharacterTable.insert(pair <string,bool> ("Cs",true));
  atomicCharacterTable.insert(pair <string,bool> ("Ba",true));
  atomicCharacterTable.insert(pair <string,bool> ("La",true));
  atomicCharacterTable.insert(pair <string,bool> ("Ce",true));
  atomicCharacterTable.insert(pair <string,bool> ("Pr",true));
  atomicCharacterTable.insert(pair <string,bool> ("Nd",true));
  atomicCharacterTable.insert(pair <string,bool> ("Pm",true));
  atomicCharacterTable.insert(pair <string,bool> ("Sm",true));
  atomicCharacterTable.insert(pair <string,bool> ("Eu",true));
  atomicCharacterTable.insert(pair <string,bool> ("Gd",true));
  atomicCharacterTable.insert(pair <string,bool> ("Tb",true));
  atomicCharacterTable.insert(pair <string,bool> ("Dy",true));
  atomicCharacterTable.insert(pair <string,bool> ("Ho",true));
  atomicCharacterTable.insert(pair <string,bool> ("Er",true));
  atomicCharacterTable.insert(pair <string,bool> ("Tm",true));
  atomicCharacterTable.insert(pair <string,bool> ("Yb",true));
  atomicCharacterTable.insert(pair <string,bool> ("Lu",true));
  atomicCharacterTable.insert(pair <string,bool> ("Hf",true));
  atomicCharacterTable.insert(pair <string,bool> ("Ta",true));
  atomicCharacterTable.insert(pair <string,bool> ("W",true));
  atomicCharacterTable.insert(pair <string,bool> ("Re",true));
  atomicCharacterTable.insert(pair <string,bool> ("Os",true));
  atomicCharacterTable.insert(pair <string,bool> ("Ir",true));
  atomicCharacterTable.insert(pair <string,bool> ("Pt",true));
  atomicCharacterTable.insert(pair <string,bool> ("Au",true));
  atomicCharacterTable.insert(pair <string,bool> ("Hg",true));
  atomicCharacterTable.insert(pair <string,bool> ("Tl",true));
  atomicCharacterTable.insert(pair <string,bool> ("Pb",true));
  atomicCharacterTable.insert(pair <string,bool> ("Bi",true));
  atomicCharacterTable.insert(pair <string,bool> ("Po",true));
  atomicCharacterTable.insert(pair <string,bool> ("At",false));
  atomicCharacterTable.insert(pair <string,bool> ("Rn",false));
  atomicCharacterTable.insert(pair <string,bool> ("Fr",true));
  atomicCharacterTable.insert(pair <string,bool> ("Ra",true));
  atomicCharacterTable.insert(pair <string,bool> ("Ac",true));
  atomicCharacterTable.insert(pair <string,bool> ("Th",true));
  atomicCharacterTable.insert(pair <string,bool> ("Pa",true));
  atomicCharacterTable.insert(pair <string,bool> ("U",true));
  atomicCharacterTable.insert(pair <string,bool> ("Np",true));
  atomicCharacterTable.insert(pair <string,bool> ("Pu",true));
  atomicCharacterTable.insert(pair <string,bool> ("Am",true));
  atomicCharacterTable.insert(pair <string,bool> ("Cm",true));
  atomicCharacterTable.insert(pair <string,bool> ("Bk",true));
  atomicCharacterTable.insert(pair <string,bool> ("Cf",true));
  atomicCharacterTable.insert(pair <string,bool> ("Es",true));
  atomicCharacterTable.insert(pair <string,bool> ("Fm",true));
  atomicCharacterTable.insert(pair <string,bool> ("Md",true));
  atomicCharacterTable.insert(pair <string,bool> ("No",true));
  atomicCharacterTable.insert(pair <string,bool> ("Lr",true));
  atomicCharacterTable.insert(pair <string,bool> ("Rf",true));
  atomicCharacterTable.insert(pair <string,bool> ("Db",true));
  atomicCharacterTable.insert(pair <string,bool> ("Sg",true));
  atomicCharacterTable.insert(pair <string,bool> ("Bh",true));
  atomicCharacterTable.insert(pair <string,bool> ("Hs",true));
  atomicCharacterTable.insert(pair <string,bool> ("Mt",true));
  atomicCharacterTable.insert(pair <string,bool> ("Ds",true));
}


/** Fills the atomic number table with atomic number of all elements
 ** */
void initializeAtomicNumberTable(){
  atomicNumberTable.insert(pair <string,int> ("H", 1));
  atomicNumberTable.insert(pair <string,int> ("D", 1));
  atomicNumberTable.insert(pair <string,int> ("He", 2));
  atomicNumberTable.insert(pair <string,int> ("Li", 3));
  atomicNumberTable.insert(pair <string,int> ("Be", 4));
  atomicNumberTable.insert(pair <string,int> ("B", 5));
  atomicNumberTable.insert(pair <string,int> ("C", 6));
  atomicNumberTable.insert(pair <string,int> ("N", 7));
  atomicNumberTable.insert(pair <string,int> ("O", 8));
  atomicNumberTable.insert(pair <string,int> ("F", 9));
  atomicNumberTable.insert(pair <string,int> ("Ne", 10));
  atomicNumberTable.insert(pair <string,int> ("Na", 11));
  atomicNumberTable.insert(pair <string,int> ("Mg", 12));
  atomicNumberTable.insert(pair <string,int> ("Al", 13));
  atomicNumberTable.insert(pair <string,int> ("Si", 14));
  atomicNumberTable.insert(pair <string,int> ("P", 15));
  atomicNumberTable.insert(pair <string,int> ("S", 16));
  atomicNumberTable.insert(pair <string,int> ("Cl", 17));
  atomicNumberTable.insert(pair <string,int> ("Ar", 18));
  atomicNumberTable.insert(pair <string,int> ("K", 19));
  atomicNumberTable.insert(pair <string,int> ("Ca", 20));
  atomicNumberTable.insert(pair <string,int> ("Sc", 21));
  atomicNumberTable.insert(pair <string,int> ("Ti", 22));
  atomicNumberTable.insert(pair <string,int> ("V", 23));
  atomicNumberTable.insert(pair <string,int> ("Cr", 24));
  atomicNumberTable.insert(pair <string,int> ("Mn", 25));
  atomicNumberTable.insert(pair <string,int> ("Fe", 26));
  atomicNumberTable.insert(pair <string,int> ("Co", 27));
  atomicNumberTable.insert(pair <string,int> ("Ni", 28));
  atomicNumberTable.insert(pair <string,int> ("Cu", 29));
  atomicNumberTable.insert(pair <string,int> ("Zn", 30));
  atomicNumberTable.insert(pair <string,int> ("Ga", 31));
  atomicNumberTable.insert(pair <string,int> ("Ge", 32));
  atomicNumberTable.insert(pair <string,int> ("As", 33));
  atomicNumberTable.insert(pair <string,int> ("Se", 34));
  atomicNumberTable.insert(pair <string,int> ("Br", 35));
  atomicNumberTable.insert(pair <string,int> ("Kr", 36));
  atomicNumberTable.insert(pair <string,int> ("Rb", 37));
  atomicNumberTable.insert(pair <string,int> ("Sr", 38));
  atomicNumberTable.insert(pair <string,int> ("Y", 39));
  atomicNumberTable.insert(pair <string,int> ("Zr", 40));
  atomicNumberTable.insert(pair <string,int> ("Nb", 41));
  atomicNumberTable.insert(pair <string,int> ("Mo", 42));
  atomicNumberTable.insert(pair <string,int> ("Tc", 43));
  atomicNumberTable.insert(pair <string,int> ("Ru", 44));
  atomicNumberTable.insert(pair <string,int> ("Rh", 45));
  atomicNumberTable.insert(pair <string,int> ("Pd", 46));
  atomicNumberTable.insert(pair <string,int> ("Ag", 47));
  atomicNumberTable.insert(pair <string,int> ("Cd", 48));
  atomicNumberTable.insert(pair <string,int> ("In", 49));
  atomicNumberTable.insert(pair <string,int> ("Sn", 50));
  atomicNumberTable.insert(pair <string,int> ("Sb", 51));
  atomicNumberTable.insert(pair <string,int> ("Te", 52));
  atomicNumberTable.insert(pair <string,int> ("I", 53));
  atomicNumberTable.insert(pair <string,int> ("Xe", 54));
  atomicNumberTable.insert(pair <string,int> ("Cs", 55));
  atomicNumberTable.insert(pair <string,int> ("Ba", 56));
  atomicNumberTable.insert(pair <string,int> ("La", 57));
  atomicNumberTable.insert(pair <string,int> ("Ce", 58));
  atomicNumberTable.insert(pair <string,int> ("Pr", 59));
  atomicNumberTable.insert(pair <string,int> ("Nd", 60));
  atomicNumberTable.insert(pair <string,int> ("Pm", 61));
  atomicNumberTable.insert(pair <string,int> ("Sm", 62));
  atomicNumberTable.insert(pair <string,int> ("Eu", 63));
  atomicNumberTable.insert(pair <string,int> ("Gd", 64));
  atomicNumberTable.insert(pair <string,int> ("Tb", 65));
  atomicNumberTable.insert(pair <string,int> ("Dy", 66));
  atomicNumberTable.insert(pair <string,int> ("Ho", 67));
  atomicNumberTable.insert(pair <string,int> ("Er", 68));
  atomicNumberTable.insert(pair <string,int> ("Tm", 69));
  atomicNumberTable.insert(pair <string,int> ("Yb", 70));
  atomicNumberTable.insert(pair <string,int> ("Lu", 71));
  atomicNumberTable.insert(pair <string,int> ("Hf", 72));
  atomicNumberTable.insert(pair <string,int> ("Ta", 73));
  atomicNumberTable.insert(pair <string,int> ("W", 74));
  atomicNumberTable.insert(pair <string,int> ("Re", 75));
  atomicNumberTable.insert(pair <string,int> ("Os", 76));
  atomicNumberTable.insert(pair <string,int> ("Ir", 77));
  atomicNumberTable.insert(pair <string,int> ("Pt", 78));
  atomicNumberTable.insert(pair <string,int> ("Au", 79));
  atomicNumberTable.insert(pair <string,int> ("Hg", 80));
  atomicNumberTable.insert(pair <string,int> ("Tl", 81));
  atomicNumberTable.insert(pair <string,int> ("Pb", 82));
  atomicNumberTable.insert(pair <string,int> ("Bi", 83));
  atomicNumberTable.insert(pair <string,int> ("Po", 84));
  atomicNumberTable.insert(pair <string,int> ("At", 85));
  atomicNumberTable.insert(pair <string,int> ("Rn", 86));
  atomicNumberTable.insert(pair <string,int> ("Fr", 87));
  atomicNumberTable.insert(pair <string,int> ("Ra", 88));
  atomicNumberTable.insert(pair <string,int> ("Ac", 89));
  atomicNumberTable.insert(pair <string,int> ("Th", 90));
  atomicNumberTable.insert(pair <string,int> ("Pa", 91));
  atomicNumberTable.insert(pair <string,int> ("U", 92));
  atomicNumberTable.insert(pair <string,int> ("Np", 93));
  atomicNumberTable.insert(pair <string,int> ("Pu", 94));
  atomicNumberTable.insert(pair <string,int> ("Am", 95));
  atomicNumberTable.insert(pair <string,int> ("Cm", 96));
  atomicNumberTable.insert(pair <string,int> ("Bk", 97));
  atomicNumberTable.insert(pair <string,int> ("Cf", 98));
  atomicNumberTable.insert(pair <string,int> ("Es", 99));
  atomicNumberTable.insert(pair <string,int> ("Fm", 100));
  atomicNumberTable.insert(pair <string,int> ("Md", 101));
  atomicNumberTable.insert(pair <string,int> ("No", 102));
  atomicNumberTable.insert(pair <string,int> ("Lr", 103));
  atomicNumberTable.insert(pair <string,int> ("Rf", 104));
  atomicNumberTable.insert(pair <string,int> ("Db", 105));
  atomicNumberTable.insert(pair <string,int> ("Sg", 106));
  atomicNumberTable.insert(pair <string,int> ("Bh", 107));
  atomicNumberTable.insert(pair <string,int> ("Hs", 108));
  atomicNumberTable.insert(pair <string,int> ("Mt", 109));
  atomicNumberTable.insert(pair <string,int> ("Ds", 110));
  atomicNumberTable.insert(pair <string,int> ("Rg", 111));
  atomicNumberTable.insert(pair <string,int> ("Cn", 112));
  atomicNumberTable.insert(pair <string,int> ("Uut", 113));
  atomicNumberTable.insert(pair <string,int> ("Fl", 114));
  atomicNumberTable.insert(pair <string,int> ("Uup", 115));
  atomicNumberTable.insert(pair <string,int> ("Lv", 116));
  atomicNumberTable.insert(pair <string,int> ("Uus", 117));
  atomicNumberTable.insert(pair <string,int> ("Uuo", 118));
}


/** Initialize Periodi Table with atom names to be used to clean atom names from indexes
 * */
void initializePT() {
        periodicTable.insert("H");
        periodicTable.insert("D");
        periodicTable.insert("He");
        periodicTable.insert("Li");
        periodicTable.insert("Be");
        periodicTable.insert("B");
        periodicTable.insert("C");
        periodicTable.insert("N");
        periodicTable.insert("O");
        periodicTable.insert("F");
        periodicTable.insert("Ne");
        periodicTable.insert("Na");
        periodicTable.insert("Mg");
        periodicTable.insert("Al");
        periodicTable.insert("Si");
        periodicTable.insert("P");
        periodicTable.insert("S");
        periodicTable.insert("Cl");
        periodicTable.insert("Ar");
        periodicTable.insert("K");
        periodicTable.insert("Ca");
        periodicTable.insert("Sc");
        periodicTable.insert("Ti");
        periodicTable.insert("V");
        periodicTable.insert("Cr");
        periodicTable.insert("Mn");
        periodicTable.insert("Fe");
        periodicTable.insert("Co");
        periodicTable.insert("Ni");
        periodicTable.insert("Cu");
        periodicTable.insert("Zn");
        periodicTable.insert("Ga");
        periodicTable.insert("Ge");
        periodicTable.insert("As");
        periodicTable.insert("Se");
        periodicTable.insert("Br");
        periodicTable.insert("Kr");
        periodicTable.insert("Rb");
        periodicTable.insert("Sr");
        periodicTable.insert("Y");
        periodicTable.insert("Zr");
        periodicTable.insert("Nb");
        periodicTable.insert("Mo");
        periodicTable.insert("Tc");
        periodicTable.insert("Ru");
        periodicTable.insert("Rh");
        periodicTable.insert("Pd");
        periodicTable.insert("Ag");
        periodicTable.insert("Cd");
        periodicTable.insert("In");
        periodicTable.insert("Sn");
        periodicTable.insert("Sb");
        periodicTable.insert("Te");
        periodicTable.insert("I");
        periodicTable.insert("Xe");
        periodicTable.insert("Cs");
        periodicTable.insert("Ba");
        periodicTable.insert("La");
        periodicTable.insert("Ce");
        periodicTable.insert("Pr");
        periodicTable.insert("Nd");
        periodicTable.insert("Pm");
        periodicTable.insert("Sm");
        periodicTable.insert("Eu");
        periodicTable.insert("Gd");
        periodicTable.insert("Tb");
        periodicTable.insert("Dy");
        periodicTable.insert("Ho");
        periodicTable.insert("Er");
        periodicTable.insert("Tm");
        periodicTable.insert("Yb");
        periodicTable.insert("Lu");
        periodicTable.insert("Hf");
        periodicTable.insert("Ta");
        periodicTable.insert("W");
        periodicTable.insert("Re");
        periodicTable.insert("Os");
        periodicTable.insert("Ir");
        periodicTable.insert("Pt");
        periodicTable.insert("Au");
        periodicTable.insert("Hg");
        periodicTable.insert("Tl");
        periodicTable.insert("Pb");
        periodicTable.insert("Bi");
        periodicTable.insert("Po");
        periodicTable.insert("At");
        periodicTable.insert("Rn");
        periodicTable.insert("Fr");
        periodicTable.insert("Ra");
        periodicTable.insert("Ac");
        periodicTable.insert("Th");
        periodicTable.insert("Pa");
        periodicTable.insert("U");
        periodicTable.insert("Np");
        periodicTable.insert("Pu");
        periodicTable.insert("Am");
        periodicTable.insert("Cm");
        periodicTable.insert("Bk");
        periodicTable.insert("Cf");
        periodicTable.insert("Es");
        periodicTable.insert("Fm");
        periodicTable.insert("Md");
        periodicTable.insert("No");
        periodicTable.insert("Lr");
        periodicTable.insert("Rf");
        periodicTable.insert("Db");
        periodicTable.insert("Sg");
        periodicTable.insert("Bh");
        periodicTable.insert("Hs");
        periodicTable.insert("Mt");
        periodicTable.insert("Ds");
}

/** Reads the radius table from the provided filename. The file must be
    formatted in two columns: atom name and radius. */
void readRadTable(char *filename){
  radTable.clear();
  fstream input;
  input.open(filename, fstream::in);
  
  if(!input.is_open()){
    cerr << "Failed to open radius input file " << filename << "\n";
    cerr << "Exiting ..." << "\n";
    exit(1);
  }
  else{
    string type = "N/A";
    double radius = -1;

    while(!input.eof()){
      input >> type >> radius;
      radTable.insert(pair<string,double> (type,radius));
    }
  }
  input.close();
}


/** Reads the mass table from the provided filename. The file must be
    formatted in two columns: atom name and mass in g/mole. */
void readMassTable(char *filename){
  massTable.clear();
  fstream input;
  input.open(filename);
  
  if(!input.is_open()){
    cerr << "Failed to open molar mass input file " << filename << "\n";
    cerr << "Exiting ..." << "\n";
    exit(1);
  }
  else{
    string type = "N/A";
    double mass = -1;

    while(!input.eof()){
      input >> type >> mass;
      massTable.insert(pair<string,double> (type,mass));
    }
  }
  input.close();
}

/** Return the radius for the corresponding atom name. If the -nor
    option was specified, returns 0. */
double lookupRadius(string atomType, bool radial){

  if(stripAtomNameInternalFlag == true) atomType = stripAtomName(atomType);

  if(!radial)
    return 0.0;
  map <string,double>::iterator  info = radTable.find(atomType);
  if(info == radTable.end()){
    cerr << "Unable to find radius for " << atomType << " in table. Please provide it " << "\n"
	 << "in a reference file or check you input file." << "\n"
	 << "Exiting ..." << "\n";
    exit(1);
  }
  else
    return info->second;
}

/** Return the covalent radius for the corresponding atom name. */
double lookupCovRadius(string atomType){
  map <string,double>::iterator  info = covRadTable.find(atomType);
  if(info == covRadTable.end()){
    cerr << "Unable to find covalent radius for " << atomType << " in table. Please modify networkinfo.cc and recomplie the code " << "\n"
         << "Exiting ..." << "\n";
    exit(1);
  }
  else
   return info->second;
}


/** Return the mass for the corresponding atom name. */
double lookupMass(string atomType){

  if(stripAtomNameInternalFlag == true) atomType = stripAtomName(atomType);

  map <string,double>::iterator  info = massTable.find(atomType);
  if(info == massTable.end()){
    cerr << "Unable to find molar mass for " << atomType << " in table. Please provide it " << "\n"
	 << "in a reference file or check you input file." << "\n"
	 << "Exiting ..." << "\n";
    exit(1);
  }
  else
    return info->second;
}


/** Return the atomic number for the corresponding atom name. */
int lookupAtomicNumber(string atomType){
  map <string,int>::iterator  info = atomicNumberTable.find(atomType);
  if(info == atomicNumberTable.end()){
    cerr << "Unable to find atomic number for " << atomType << " in table. Please provide it " << "\n"
         << "in the source code and recompile the code." << "\n"
         << "Exiting ..." << "\n";
    exit(1);
  }
  else
    return info->second;
}

/** Return true is an atom is a metal **/
bool isMetal(string atomType){
  map <string,bool>::iterator  info = atomicCharacterTable.find(atomType);
  if(info == atomicCharacterTable.end()){
    cerr << "Unable to find character information for " << atomType << " in table. Please modify networkinfo.cc and recomplie the code " << "\n" 
         << "Exiting ..." << "\n";
    exit(1);
  }
  else
    return info->second;
}

/** Clean up atom name to remove any index strings */
std::string stripAtomName(std::string extAtom) {
        string atom2 = extAtom.substr(0, 2);
        string atom1 = extAtom.substr(0, 1);
        if (periodicTable.find(atom2) != periodicTable.end()) {
                return atom2;
        } else if (periodicTable.find(atom1) != periodicTable.end()) {
                return atom1;
        } else {
                return extAtom;
        }
}

